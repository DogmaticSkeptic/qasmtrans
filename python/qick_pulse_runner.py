#!/usr/bin/env python3
"""Utility for sending QASMTrans pulse schedules to a QICK board.

This script loads the pulse schedule JSON produced by QASMTrans's
``dumpPulses`` routine, maps each pulse to a QICK signal generator, and builds
an ``AveragerProgram`` that replays the sequence on hardware.  It is intended
as a starting point: the mapping from logical qubits to physical generator and
readout channels is provided via a simple configuration structure that you can
adjust for your lab setup.

Example usage (dry run):

    python qick_pulse_runner.py data/output/my_transpiled_pulses.json \
        --config qick_config.json --summary

To execute on hardware (requires qick, pynq, etc.):

    python qick_pulse_runner.py data/output/my_transpiled_pulses.json \
        --config qick_config.json --run --decimated
"""
from __future__ import annotations

import argparse
import json
import logging
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional

try:  # Guard import so the script can still print summaries without qick.
    from qick import AveragerProgram, QickConfig, QickSoc
    QICK_AVAILABLE = True
except ImportError:  # pragma: no cover - optional dependency for planning
    AveragerProgram = None  # type: ignore
    QickConfig = object  # type: ignore[assignment]
    QickSoc = None  # type: ignore
    QICK_AVAILABLE = False

PULSE_GAIN_MAX = 2 ** 15 - 1  # QICK DAC gain register is 16-bit signed


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def normalise_qubit_map(raw_map: Dict[str, int]) -> Dict[str, int]:
    normalised: Dict[str, int] = {}
    for key, value in raw_map.items():
        try:
            qubit = str(int(key))
            normalised[qubit] = int(value)
        except (ValueError, TypeError) as exc:
            raise ValueError(f"Invalid qubit mapping entry {key!r}: {value!r}") from exc
    return normalised


def default_config() -> dict:
    return {
        "qubit_gen_map": {"0": 6},  # logical qubit -> QICK generator channel
        "ro_chs": [0],
        "pulse_freq": 250.0,  # MHz
        "readout_freq": 250.0,  # MHz (can be overridden per readout)
        "pulse_gain_scale": 30000,  # DAC units per unit amplitude
        "gain_max": PULSE_GAIN_MAX,
        "init_synci": 200,
        "relax_delay": 1.0,  # microseconds
        "adc_trig_offset": 100,
        "readout_length": 200,
        "soft_avgs": 10,
        "reps": 1,
        "enable_readout": False,
        "nqz": 1,
        "start_offset": 0.0,  # seconds
    }


def merge_user_config(base: dict, override: Optional[Path]) -> dict:
    config = dict(base)
    if override is None:
        return config
    user_cfg = load_json(override)
    for key, value in user_cfg.items():
        config[key] = value
    config["qubit_gen_map"] = normalise_qubit_map(config.get("qubit_gen_map", {}))
    config["ro_chs"] = [int(ch) for ch in config.get("ro_chs", [])]
    config.setdefault("pulse_gain_scale", 30000)
    config.setdefault("gain_max", PULSE_GAIN_MAX)
    config.setdefault("enable_readout", False)
    config.setdefault("relax_delay", 1.0)
    config.setdefault("readout_length", 200)
    config.setdefault("soft_avgs", 1)
    config.setdefault("reps", 1)
    config.setdefault("nqz", 1)
    config.setdefault("start_offset", 0.0)
    return config


@dataclass
class PulseEvent:
    gate: str
    pulse_id: str
    qubits: List[int]
    gen_ch: int
    start_time_s: float
    width_s: float
    amplitude: float
    shape: str
    parameters: dict
    note: str


def prepare_pulse_events(pulse_doc: dict, config: dict) -> List[PulseEvent]:
    library = {entry["id"]: entry for entry in pulse_doc.get("pulse_library", [])}
    schedule = pulse_doc.get("schedule", [])
    events: List[PulseEvent] = []
    skipped_without_mapping: set[str] = set()
    skipped_zero_width = 0
    for item in schedule:
        pulse_id = item.get("pulse_id")
        definition = library.get(pulse_id)
        if definition is None:
            logging.warning("Pulse %s missing from library; skipping", pulse_id)
            continue
        width = float(definition.get("width", 0.0))
        if width <= 0.0:
            skipped_zero_width += 1
            continue
        qubits = item.get("qubits", [])
        if not qubits:
            logging.debug("Skipping gate %s with no qubits", pulse_id)
            continue
        driver_qubit = str(qubits[0])
        if driver_qubit not in config["qubit_gen_map"]:
            skipped_without_mapping.add(driver_qubit)
            continue
        gen_ch = config["qubit_gen_map"][driver_qubit]
        parameters = item.get("parameters", {})
        events.append(
            PulseEvent(
                gate=item.get("gate", ""),
                pulse_id=pulse_id,
                qubits=list(qubits),
                gen_ch=int(gen_ch),
                start_time_s=float(item.get("start_time", 0.0)),
                width_s=width,
                amplitude=float(definition.get("amplitude", 0.0)),
                shape=str(definition.get("shape", "")),
                parameters=dict(parameters),
                note=str(definition.get("note", "")),
            )
        )
    if skipped_without_mapping:
        logging.warning(
            "Skipped %d pulses – no generator mapping for logical qubits: %s",
            len(skipped_without_mapping),
            ", ".join(sorted(skipped_without_mapping)),
        )
    if skipped_zero_width:
        logging.info("Ignored %d zero-width pulses (virtual rotations)", skipped_zero_width)
    events.sort(key=lambda event: event.start_time_s)
    config["generators"] = sorted({event.gen_ch for event in events})
    return events


if QICK_AVAILABLE:

    class TranspiledPulseProgram(AveragerProgram):
        """Simple AveragerProgram that replays a fixed pulse schedule."""

        def __init__(self, soccfg: QickConfig, cfg: dict, events: List[PulseEvent]):
            self.events = events
            super().__init__(soccfg, cfg)

        def initialize(self):
            cfg = self.cfg
            self.compiled_events: List[dict] = []
            generators = cfg.get("generators", [])
            if not generators:
                raise ValueError("No pulse events to schedule. Nothing to do.")
            for gen_ch in generators:
                self.declare_gen(ch=gen_ch, nqz=cfg.get("nqz", 1))
            readout_channels = cfg.get("ro_chs", [])
            readout_freq = cfg.get("readout_freq", cfg.get("pulse_freq"))
            readout_gen = cfg.get("readout_gen_ch", generators[0])
            if cfg.get("enable_readout", False):
                if readout_freq is None:
                    raise ValueError("Readout requested but readout_freq is None")
                for ro_ch in readout_channels:
                    self.declare_readout(
                        ch=ro_ch,
                        length=int(cfg.get("readout_length", 200)),
                        freq=readout_freq,
                        gen_ch=readout_gen,
                    )
            default_phase = cfg.get("default_phase_deg", 0.0)
            qubit_freqs = {str(k): v for k, v in cfg.get("qubit_freqs", {}).items()}
            for qubit, gen_ch in cfg.get("qubit_gen_map", {}).items():
                freq = qubit_freqs.get(str(qubit), cfg.get("pulse_freq"))
                if freq is None:
                    continue
                freq_reg = self.freq2reg(
                    freq, gen_ch=gen_ch, ro_ch=readout_channels[0] if readout_channels else None
                )
                phase_reg = self.deg2reg(default_phase, gen_ch=gen_ch)
                default_gain = int(cfg.get("default_gain", cfg.get("pulse_gain_scale", 30000)))
                default_gain = max(0, min(cfg.get("gain_max", PULSE_GAIN_MAX), default_gain))
                self.default_pulse_registers(ch=gen_ch, freq=freq_reg, phase=phase_reg, gain=default_gain)
            gain_scale = cfg.get("pulse_gain_scale", 30000)
            gain_max = cfg.get("gain_max", PULSE_GAIN_MAX)
            start_offset = float(cfg.get("start_offset", 0.0))
            for event in self.events:
                length_cycles = max(
                    1, int(round(self.us2cycles(event.width_s * 1e6, gen_ch=event.gen_ch)))
                )
                start_cycles = int(round(self.us2cycles((event.start_time_s + start_offset) * 1e6)))
                gain = int(round(event.amplitude * gain_scale))
                gain = max(0, min(gain_max, gain))
                compiled = {
                    "gate": event.gate,
                    "pulse_id": event.pulse_id,
                    "gen_ch": event.gen_ch,
                    "length": length_cycles,
                    "t_start": start_cycles,
                    "gain": gain,
                    "parameters": event.parameters,
                    "note": event.note,
                }
                self.compiled_events.append(compiled)
            self.compiled_events.sort(key=lambda entry: entry["t_start"])
            self.synci(cfg.get("init_synci", 200))

        def body(self):
            cfg = self.cfg
            if not self.compiled_events:
                self.sync_all(self.us2cycles(cfg.get("relax_delay", 1.0)))
                return
            for entry in self.compiled_events:
                if entry["length"] <= 0:
                    continue
                self.set_pulse_registers(
                    ch=entry["gen_ch"],
                    style="const",
                    length=entry["length"],
                    gain=entry["gain"],
                )
                self.pulse(ch=entry["gen_ch"], t=entry["t_start"])
            self.wait_all()
            self.sync_all(self.us2cycles(cfg.get("relax_delay", 1.0)))

else:  # pragma: no cover - executed only when qick is absent

    class TranspiledPulseProgram:  # type: ignore[override]
        def __init__(self, *args, **kwargs):
            raise RuntimeError(
                "qick package is not available. Install qick or rerun with --summary to avoid execution."
            )


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__ or "")
    parser.add_argument("pulse_file", type=Path, help="Pulse JSON emitted by QASMTrans")
    parser.add_argument(
        "--config", type=Path, default=None, help="Optional JSON config mapping qubits to QICK channels"
    )
    parser.add_argument(
        "--summary",
        action="store_true",
        help="Print a human readable summary of the loaded pulses and exit",
    )
    parser.add_argument(
        "--run",
        action="store_true",
        help="Execute the sequence on connected QICK hardware",
    )
    parser.add_argument(
        "--decimated",
        action="store_true",
        help="Use acquire_decimated() when running (default is acquire())",
    )
    parser.add_argument(
        "--progress",
        action="store_true",
        help="Show acquisition progress bar when running with acquire_decimated",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"],
        help="Logging verbosity",
    )
    return parser


def summarise_events(events: List[PulseEvent]) -> str:
    if not events:
        return "No eligible pulses found in the schedule."
    lines = [f"Prepared {len(events)} pulses mapped to generators: {sorted({e.gen_ch for e in events})}"]
    by_gate: Dict[str, int] = {}
    for event in events:
        by_gate[event.gate] = by_gate.get(event.gate, 0) + 1
    gate_summary = ", ".join(f"{gate}:{count}" for gate, count in sorted(by_gate.items()))
    lines.append(f"Gate histogram: {gate_summary}")
    total_time = events[-1].start_time_s + events[-1].width_s
    lines.append(f"Total scheduled duration: {total_time * 1e6:.3f} us")
    return "\n".join(lines)


def ensure_qick_available():
    if AveragerProgram is None or QickSoc is None:
        raise RuntimeError(
            "qick package is not available. Install qick or run with --summary for an offline preview."
        )


def run_program(pulse_file: Path, config_path: Optional[Path], run: bool, decimated: bool, progress: bool):
    pulse_doc = load_json(pulse_file)
    config = merge_user_config(default_config(), config_path)
    events = prepare_pulse_events(pulse_doc, config)
    print(summarise_events(events))
    if not run:
        return
    ensure_qick_available()
    if not events:
        print("No pulses to execute; exiting.")
        return
    soc = QickSoc()
    soccfg = soc  # local execution; adjust if running remotely
    program = TranspiledPulseProgram(soccfg, config, events)
    if decimated:
        iq_data = program.acquire_decimated(soc, progress=progress)
        print(f"Collected {len(iq_data)} decimated traces with shape {[len(trace[0]) for trace in iq_data]} samples")
    else:
        avgi, avgq = program.acquire(soc)
        print(f"Accumulated I averages: {avgi}")
        print(f"Accumulated Q averages: {avgq}")


def main(argv: Optional[Iterable[str]] = None) -> int:
    parser = build_argument_parser()
    args = parser.parse_args(list(argv) if argv is not None else None)
    logging.basicConfig(level=getattr(logging, args.log_level))
    try:
        run_program(
            pulse_file=args.pulse_file,
            config_path=args.config,
            run=args.run,
            decimated=args.decimated,
            progress=args.progress,
        )
    except Exception as exc:  # noqa: BLE001 - present concise error to CLI
        logging.error("%s", exc)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
