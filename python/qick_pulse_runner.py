#!/usr/bin/env python3
"""Play pulse schedules produced by QASMTrans on a QICK board."""
import argparse
import json
import sys
from collections import Counter
from pathlib import Path

from qick import AveragerProgram, QickSoc

PULSE_GAIN_MAX = (1 << 15) - 1


def load_json(path):
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def normalise_qubit_map(raw_map):
    return {str(int(k)): int(v) for k, v in raw_map.items()}


def merge_user_config(path):
    cfg = {
        "qubit_gen_map": {"0": 6},
        "ro_chs": [0],
        "pulse_freq": 250.0,
        "readout_freq": 250.0,
        "pulse_gain_scale": 30000,
        "gain_max": PULSE_GAIN_MAX,
        "init_synci": 200,
        "relax_delay": 1.0,
        "adc_trig_offset": 100,
        "readout_length": 200,
        "soft_avgs": 10,
        "reps": 1,
        "enable_readout": False,
        "nqz": 1,
        "start_offset": 0.0,
        "scope_pin": True,
    }
    if path is not None:
        cfg.update(load_json(path))
    cfg["qubit_gen_map"] = normalise_qubit_map(cfg.get("qubit_gen_map", {}))
    cfg["ro_chs"] = [int(ch) for ch in cfg.get("ro_chs", [])]
    return cfg


def prepare_events(pulse_doc, cfg):
    lib = {p["id"]: p for p in pulse_doc.get("pulse_library", [])}
    events = []
    for item in pulse_doc.get("schedule", []):
        qubits = item.get("qubits") or []
        pulse = lib.get(item.get("pulse_id"))
        if not qubits or pulse is None or float(pulse.get("width", 0)) <= 0:
            continue
        driver = str(qubits[0])
        if driver not in cfg["qubit_gen_map"]:
            continue
        events.append(
            {
                "gate": item.get("gate", ""),
                "pulse_id": item.get("pulse_id"),
                "qubits": [int(q) for q in qubits],
                "gen_ch": int(cfg["qubit_gen_map"][driver]),
                "start_time_s": float(item.get("start_time", 0.0)),
                "width_s": float(pulse.get("width", 0.0)),
                "amplitude": float(pulse.get("amplitude", 0.0)),
            }
        )
    events.sort(key=lambda e: e["start_time_s"])
    cfg["generators"] = sorted({e["gen_ch"] for e in events})
    return events


class TranspiledPulseProgram(AveragerProgram):
    def __init__(self, soccfg, cfg, events):
        self.events = events
        super().__init__(soccfg, cfg)

    def initialize(self):
        cfg = self.cfg
        for ch in cfg["generators"]:
            self.declare_gen(ch=ch, nqz=cfg.get("nqz", 1))
        if cfg.get("enable_readout") and cfg["ro_chs"]:
            readout_gen = cfg.get("readout_gen_ch", cfg["generators"][0])
            for ro in cfg["ro_chs"]:
                self.declare_readout(
                    ch=ro,
                    length=int(cfg.get("readout_length", 200)),
                    freq=cfg.get("readout_freq", cfg.get("pulse_freq")),
                    gen_ch=readout_gen,
                )
        self.specs = []
        gain_scale = cfg.get("pulse_gain_scale", 30000)
        gain_max = cfg.get("gain_max", PULSE_GAIN_MAX)
        start_offset = float(cfg.get("start_offset", 0.0))
        default_phase = cfg.get("default_phase_deg", 0.0)
        qubit_freqs = {str(k): v for k, v in cfg.get("qubit_freqs", {}).items()}
        for qubit, ch in cfg.get("qubit_gen_map", {}).items():
            freq = qubit_freqs.get(qubit, cfg.get("pulse_freq"))
            if freq is None:
                continue
            self.default_pulse_registers(
                ch=ch,
                freq=self.freq2reg(freq, gen_ch=ch, ro_ch=cfg["ro_chs"][0] if cfg["ro_chs"] else None),
                phase=self.deg2reg(default_phase, gen_ch=ch),
                gain=min(gain_max, int(cfg.get("default_gain", gain_scale))),
            )
        for event in self.events:
            length = max(1, int(round(self.us2cycles(event["width_s"] * 1e6, gen_ch=event["gen_ch"]))))
            start = int(round(self.us2cycles((event["start_time_s"] + start_offset) * 1e6)))
            gain = min(gain_max, int(round(event["amplitude"] * gain_scale)))
            self.specs.append({"gen_ch": event["gen_ch"], "start": start, "length": length, "gain": gain})
        self.specs.sort(key=lambda spec: spec["start"])
        self.synci(cfg.get("init_synci", 200))

    def body(self):
        cfg = self.cfg
        if cfg.get("enable_readout") and cfg["ro_chs"]:
            self.trigger(
                adcs=cfg["ro_chs"],
                pins=[0] if cfg.get("scope_pin", True) else [],
                adc_trig_offset=cfg.get("adc_trig_offset", 100),
            )
        for spec in self.specs:
            self.set_pulse_registers(ch=spec["gen_ch"], style="const", length=spec["length"], gain=spec["gain"])
            self.pulse(ch=spec["gen_ch"], t=spec["start"])
        self.wait_all()
        self.sync_all(self.us2cycles(cfg.get("relax_delay", 1.0)))


def run_program(pulse_file, config_path, run, decimated, progress):
    cfg = merge_user_config(config_path)
    events = prepare_events(load_json(pulse_file), cfg)
    if events:
        counts = Counter(event["gate"] for event in events)
        total = events[-1]["start_time_s"] + events[-1]["width_s"]
        print(
            f"Prepared {len(events)} pulses mapped to generators: {sorted({e['gen_ch'] for e in events})}\n"
            f"Gate histogram: {', '.join(f'{g}:{c}' for g, c in sorted(counts.items()))}\n"
            f"Total scheduled duration: {total * 1e6:.3f} us"
        )
    else:
        print("No eligible pulses found in the schedule.")
    if not run or not events:
        return
    soc = QickSoc()
    program = TranspiledPulseProgram(soc, cfg, events)
    if decimated:
        traces = program.acquire_decimated(soc, progress=progress)
        print("Collected decimated traces: " + ", ".join(str(len(t[0])) for t in traces))
    else:
        avgi, avgq = program.acquire(soc)
        print(f"Accumulated I averages: {avgi}")
        print(f"Accumulated Q averages: {avgq}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__ or "")
    parser.add_argument("pulse_file", type=Path)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--summary", action="store_true")
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--decimated", action="store_true")
    parser.add_argument("--progress", action="store_true")
    args = parser.parse_args(list(argv) if argv is not None else None)
    run_program(args.pulse_file, args.config, args.run and not args.summary, args.decimated, args.progress)


if __name__ == "__main__":
    sys.exit(main())
