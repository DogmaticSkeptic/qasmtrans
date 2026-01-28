#!/usr/bin/env python3
"""Play pulse schedules produced by QASMTrans on a QICK board."""
import argparse
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

try:
    from qick import AveragerProgram, QickSoc
except ImportError:  # pragma: no cover - optional dependency for summary-only usage
    AveragerProgram = None
    QickSoc = None

PULSE_GAIN_MAX = (1 << 15) - 1


def load_json(path):
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _prepare_waveform(samples_i, samples_q, gain_max):
    if not samples_i and not samples_q:
        return [], []
    combined = list(samples_i) + list(samples_q)
    max_abs = max((abs(v) for v in combined), default=1.0) or 1.0
    scale = gain_max / max_abs
    def clamp(value):
        return int(max(min(round(value * scale), gain_max), -gain_max))

    idata = [clamp(v) for v in samples_i] if samples_i else []
    qdata = [clamp(v) for v in samples_q] if samples_q else [0] * len(idata)
    if not idata and samples_q:
        idata = [0] * len(qdata)
    if len(qdata) != len(idata):
        length = max(len(idata), len(qdata))
        idata = (idata + [0] * length)[:length]
        qdata = (qdata + [0] * length)[:length]
    return idata, qdata


def prepare_events(pulse_doc, cfg):
    lib = {p["id"]: p for p in pulse_doc.get("pulse_library", [])}
    qubit_gen_map = cfg.get("qubit_gen_map", {})
    events = []
    frame_phases = defaultdict(float)
    for item in pulse_doc.get("schedule", []):
        qubits = item.get("qubits") or []
        if not qubits:
            continue
        pulse = lib.get(item.get("pulse_id"))
        if pulse is None:
            continue
        waveform_type = str(pulse.get("waveform_type", "")).lower()
        is_virtual = bool(pulse.get("virtual")) or waveform_type == "virtual" or bool(item.get("virtual"))
        if is_virtual:
            theta = None
            parameters = item.get("parameters")
            if isinstance(parameters, dict):
                theta = parameters.get("theta")
            if theta is None:
                theta = (pulse.get("parameters") or {}).get("theta")
            if theta is not None:
                theta_value = float(theta)
                for q in qubits:
                    frame_phases[str(q)] += theta_value
            continue
        duration = item.get("duration", pulse.get("width", 0.0))
        if float(duration) <= 0:
            continue
        channels = []
        missing_channel = False
        for q in qubits:
            key = str(q)
            if key not in qubit_gen_map:
                missing_channel = True
                break
            channels.append(int(qubit_gen_map[key]))
        if missing_channel or not channels:
            continue
        samples_i = [float(v) for v in pulse.get("samples_i", [])]
        samples_q = [float(v) for v in pulse.get("samples_q", [])]
        if samples_i or samples_q:
            peak_i = max((abs(v) for v in samples_i), default=0.0)
            peak_q = max((abs(v) for v in samples_q), default=0.0)
            amplitude = max(peak_i, peak_q)
        else:
            amplitude = float(pulse.get("amplitude", 0.0))
        phase_shifts_deg = [math.degrees(frame_phases[str(q)]) for q in qubits]
        events.append(
            {
                "gate": item.get("gate", ""),
                "pulse_id": item.get("pulse_id"),
                "qubits": qubits,
                "gen_ch": channels[0],
                "gen_chs": channels,
                "start_time_s": float(item.get("start_time", 0.0)),
                "width_s": float(duration),
                "amplitude": amplitude,
                "shape": pulse.get("shape"),
                "waveform_type": pulse.get("waveform_type", pulse.get("shape")),
                "samples_i": samples_i,
                "samples_q": samples_q,
                "phase_shifts_deg": phase_shifts_deg,
            }
        )
    events.sort(key=lambda e: e["start_time_s"])
    cfg["generators"] = sorted({ch for event in events for ch in event.get("gen_chs", [event.get("gen_ch")]) if ch is not None})
    return events


if AveragerProgram is None:
    class TranspiledPulseProgram:  # type: ignore
        def __init__(self, *_, **__):
            raise ImportError("qick package is required for hardware execution but is not installed.")
else:
    class TranspiledPulseProgram(AveragerProgram):
        def __init__(self, soccfg, cfg, events):
            self.events = events
            super().__init__(soccfg, cfg)

        def initialize(self):
            cfg = self.cfg
            for ch in cfg.get("generators", []):
                self.declare_gen(ch=ch, nqz=cfg.get("nqz", 1))
            if cfg.get("enable_readout") and cfg.get("ro_chs"):
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
            channel_base_phase = {}
            for qubit, ch in cfg.get("qubit_gen_map", {}).items():
                freq = qubit_freqs.get(str(qubit), cfg.get("pulse_freq"))
                if freq is None:
                    continue
                ro_channel = cfg.get("ro_chs", [None])[0]
                self.default_pulse_registers(
                    ch=ch,
                    freq=self.freq2reg(freq, gen_ch=ch, ro_ch=ro_channel),
                    phase=self.deg2reg(default_phase, gen_ch=ch),
                    gain=min(gain_max, int(cfg.get("default_gain", gain_scale))),
                )
                channel_base_phase[int(ch)] = default_phase
            loaded_waveforms = set()
            for event in self.events:
                channel_specs = []
                channels = event.get("gen_chs") or [event.get("gen_ch")]
                waveform_type = (event.get("waveform_type") or "").lower()
                samples_i = event.get("samples_i") or []
                samples_q = event.get("samples_q") or []
                phase_shifts = event.get("phase_shifts_deg") or []
                for idx, ch in enumerate(channels):
                    if ch is None:
                        continue
                    length = max(1, int(round(self.us2cycles(event["width_s"] * 1e6, gen_ch=ch))))
                    start = int(round(self.us2cycles((event["start_time_s"] + start_offset) * 1e6, gen_ch=ch)))
                    gain = min(gain_max, int(round(event["amplitude"] * gain_scale)))
                    shift = phase_shifts[idx] if idx < len(phase_shifts) else (phase_shifts[-1] if phase_shifts else 0.0)
                    base_phase = channel_base_phase.get(int(ch), default_phase)
                    total_phase = base_phase + shift
                    # Keep the phase within 0..360 for numerical stability
                    total_phase = (total_phase + 360.0) % 360.0
                    phase_reg = self.deg2reg(total_phase, gen_ch=ch)
                    waveform_name = None
                    style = "const"
                    if waveform_type in {"arbitrary", "arb"} and samples_i:
                        waveform_name = f"{event['pulse_id']}_ch{ch}"
                        cache_key = (ch, waveform_name)
                        if cache_key not in loaded_waveforms:
                            idata, qdata = _prepare_waveform(samples_i, samples_q, gain_max)
                            self.add_pulse(ch=ch, name=waveform_name, idata=idata, qdata=qdata)
                            loaded_waveforms.add(cache_key)
                        style = "arb"
                    elif waveform_type in {"flat_top", "flat-top"} and samples_i:
                        style = "flat_top"
                        waveform_name = f"{event['pulse_id']}_ch{ch}_flat"
                        cache_key = (ch, waveform_name)
                        if cache_key not in loaded_waveforms:
                            idata, qdata = _prepare_waveform(samples_i, samples_q, gain_max)
                            self.add_pulse(ch=ch, name=waveform_name, idata=idata, qdata=qdata)
                            loaded_waveforms.add(cache_key)
                    else:
                        style = "const"
                    channel_specs.append(
                        {
                            "ch": ch,
                            "start": start,
                            "length": length,
                            "gain": gain,
                            "style": style,
                            "waveform": waveform_name,
                            "phase": phase_reg,
                        }
                    )
                if not channel_specs:
                    continue
                self.specs.append({
                    "channels": channel_specs,
                    "style": channel_specs[0]["style"] if channel_specs else "const",
                })
            self.specs.sort(key=lambda spec: min(ch_spec["start"] for ch_spec in spec["channels"]))
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
                for ch_spec in spec["channels"]:
                    kwargs = {
                        "ch": ch_spec["ch"],
                        "style": ch_spec.get("style", spec.get("style", "const")),
                        "gain": ch_spec["gain"],
                    }
                    if ch_spec.get("length") and ch_spec.get("style") != "arb":
                        kwargs["length"] = ch_spec["length"]
                    if ch_spec.get("waveform"):
                        kwargs["waveform"] = ch_spec["waveform"]
                    if ch_spec.get("phase") is not None:
                        kwargs["phase"] = ch_spec["phase"]
                    self.set_pulse_registers(**kwargs)
                for ch_spec in spec["channels"]:
                    self.pulse(ch=ch_spec["ch"], t=ch_spec["start"])
            self.wait_all()
            self.sync_all(self.us2cycles(cfg.get("relax_delay", 1.0)))


def run_program(pulse_file, qick_config_path, run_enabled, summary_only):
    pulse_doc = load_json(pulse_file)
    cfg = load_json(qick_config_path)
    events = prepare_events(pulse_doc, cfg)
    if events:
        counts = Counter(event["gate"] for event in events)
        total = events[-1]["start_time_s"] + events[-1]["width_s"]
        generator_set = sorted({ch for event in events for ch in event.get("gen_chs", [event.get("gen_ch")])})
        print(
            f"Prepared {len(events)} pulses mapped to generators: {generator_set}\n"
            f"Gate histogram: {', '.join(f'{g}:{c}' for g, c in sorted(counts.items()))}\n"
            f"Total scheduled duration: {total * 1e6:.3f} us"
        )
    else:
        print("No eligible pulses found in the schedule.")
    if summary_only or not run_enabled or not events:
        return
    if QickSoc is None or AveragerProgram is None:
        raise ImportError("qick package is required for hardware execution but is not installed.")
    try:
        soc = QickSoc()
    except Exception as exc:  # pragma: no cover - hardware dependent
        raise RuntimeError(f"Unable to initialize QICK hardware (QickSoc): {exc}. "
                           "Ensure the QICK device is connected/configured.") from exc
    program = TranspiledPulseProgram(soc, cfg, events)
    avgi, avgq = program.acquire(soc)
    print(f"Accumulated I averages: {avgi}")
    print(f"Accumulated Q averages: {avgq}")


def run_program_obj(pulse_doc, qick_config, run_enabled=False, summary_only=False):
    """In-memory variant: accepts dicts instead of file paths."""
    events = prepare_events(pulse_doc, qick_config)
    summary = {}
    if events:
        counts = Counter(event["gate"] for event in events)
        total = events[-1]["start_time_s"] + events[-1]["width_s"]
        generator_set = sorted({ch for event in events for ch in event.get("gen_chs", [event.get("gen_ch")])})
        summary = {
            "event_count": len(events),
            "generators": generator_set,
            "gate_histogram": {g: c for g, c in sorted(counts.items())},
            "total_duration_us": total * 1e6,
        }
    if summary_only or not run_enabled or not events:
        return summary
    if QickSoc is None or AveragerProgram is None:
        raise ImportError("qick package is required for hardware execution but is not installed.")
    try:
        soc = QickSoc()
    except Exception as exc:  # pragma: no cover - hardware dependent
        raise RuntimeError(f"Unable to initialize QICK hardware (QickSoc): {exc}. "
                           "Ensure the QICK device is connected/configured.") from exc
    program = TranspiledPulseProgram(soc, qick_config, events)
    avgi, avgq = program.acquire(soc)
    summary["avgi"] = avgi
    summary["avgq"] = avgq
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__ or "")
    parser.add_argument("circuit_pulses", type=Path)
    parser.add_argument("qick_config", type=Path)
    parser.add_argument("--summary", action="store_true")
    parser.add_argument("--run", action="store_true")
    args = parser.parse_args(list(argv) if argv is not None else None)
    run_program(args.circuit_pulses, args.qick_config, args.run, args.summary)


if __name__ == "__main__":
    sys.exit(main())
