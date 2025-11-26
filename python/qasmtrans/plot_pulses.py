#!/usr/bin/env python3
"""Visualize circuit pulses exported by QASMTrans as DAC waveforms."""

import argparse
import json
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import matplotlib.pyplot as plt
import numpy as np


def load_json(path: Path):
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def build_qubit_events(pulse_doc):
    library = {entry["id"]: entry for entry in pulse_doc.get("pulse_library", [])}
    analog_events_by_qubit: Dict[str, List[dict]] = defaultdict(list)
    virtual_events_by_qubit: Dict[str, List[dict]] = defaultdict(list)
    label_sequences: Dict[str, List[dict]] = defaultdict(list)
    for item in pulse_doc.get("schedule", []):
        pulse_id = item.get("pulse_id")
        pulse = library.get(pulse_id, {})
        qubits = item.get("qubits") or []
        if not qubits:
            continue
        start = float(item.get("start_time", 0.0))
        duration = float(item.get("duration", pulse.get("width", 0.0)))
        gate_name = item.get("gate") or pulse.get("gate") or pulse_id or "pulse"
        label = gate_name
        is_virtual = pulse.get("virtual") or str(pulse.get("waveform_type", "")).lower() == "virtual"
        parameters = pulse.get("parameters", {})
        theta = parameters.get("theta") or item.get("parameters", {}).get("theta")
        for qubit in qubits:
            bucket = virtual_events_by_qubit if is_virtual else analog_events_by_qubit
            bucket[str(qubit)].append(
                {
                    "start_us": start * 1e6,
                    "width_us": duration * 1e6,
                    "label": f"{label}\n({pulse_id})" if pulse_id else label,
                    "waveform_type": (pulse.get("waveform_type") or "virtual"),
                    "shape": (pulse.get("shape") or "virtual"),
                    "samples_i": [float(v) for v in pulse.get("samples_i", [])],
                    "samples_q": [float(v) for v in pulse.get("samples_q", [])],
                    "amplitude": float(pulse.get("amplitude", 0.0)),
                    "theta": theta,
                    "virtual": is_virtual,
                    "gate": gate_name.lower(),
                }
            )
            label_sequences[str(qubit)].append(
                {
                    "label": f"{label}\n({pulse_id})" if pulse_id else label,
                    "theta": theta,
                    "virtual": is_virtual,
                    "gate": gate_name.lower(),
                }
            )
    for events in analog_events_by_qubit.values():
        events.sort(key=lambda entry: entry["start_us"])
    for events in virtual_events_by_qubit.values():
        events.sort(key=lambda entry: entry["start_us"])
    return dict(analog_events_by_qubit), dict(virtual_events_by_qubit), dict(label_sequences)


def sample_waveform(event: dict, samples_per_us: float):
    width_us = max(event["width_us"], 1e-6)
    samples_i = event.get("samples_i") or []
    samples_q = event.get("samples_q") or []

    if samples_i:
        i_array = np.array(samples_i, dtype=float)
        q_array = np.array(samples_q, dtype=float) if samples_q else np.zeros_like(i_array)
        sample_count = len(i_array)
        ts = np.linspace(0.0, width_us, sample_count, endpoint=False)
        return ts, i_array, q_array

    sample_count = max(32, int(width_us * samples_per_us))
    ts = np.linspace(0.0, width_us, sample_count, endpoint=False)
    normalized_t = ts / width_us if width_us > 0 else ts
    shape_lc = (event.get("waveform_type") or event.get("shape") or "const").lower()

    if shape_lc == "gaussian":
        sigma = 0.2
        envelope = np.exp(-0.5 * ((normalized_t - 0.5) / sigma) ** 2)
    elif shape_lc in {"const", "constant", "square"}:
        envelope = np.ones_like(ts)
    elif shape_lc == "flat_top":
        envelope = np.ones_like(ts)
        edge_samples = max(1, int(0.1 * len(ts)))
        if 2 * edge_samples < len(ts):
            ramp = 0.5 * (1 - np.cos(np.linspace(0.0, np.pi, 2 * edge_samples)))
            envelope[:edge_samples] = ramp[:edge_samples]
            envelope[-edge_samples:] = ramp[edge_samples:]
    elif shape_lc in {"cos", "cosine"}:
        envelope = 0.5 * (1 - np.cos(np.pi * normalized_t))
    elif shape_lc in {"hann", "hanning"}:
        envelope = 0.5 - 0.5 * np.cos(2 * np.pi * normalized_t)
    else:
        envelope = np.ones_like(ts)

    amplitude = event["amplitude"]
    return ts, amplitude * envelope, np.zeros_like(envelope)


def format_theta(theta: float) -> str:
    if theta is None:
        return "?"
    return f"{theta:.3f}"


def format_theta_compact(theta: float) -> str:
    if theta is None:
        return "?"
    if np.isclose(theta, np.pi / 2, atol=1e-6):
        return "pi/2"
    if np.isclose(theta, -np.pi / 2, atol=1e-6):
        return "-pi/2"
    if np.isclose(theta, np.pi, atol=1e-6):
        return "pi"
    if np.isclose(theta, -np.pi, atol=1e-6):
        return "-pi"
    return f"{theta:.3f}"


def plot_events(
    analog_events_by_qubit,
    virtual_events_by_qubit,
    label_sequences,
    title: Optional[str],
    output_path: Path,
    dpi: int,
    samples_per_us: float,
):
    if not analog_events_by_qubit and not virtual_events_by_qubit:
        raise ValueError("No pulse events were found in the schedule.")

    qubit_ids = sorted(
        set(analog_events_by_qubit.keys()) | set(virtual_events_by_qubit.keys()),
        key=lambda q: int(q) if str(q).isdigit() else str(q),
    )
    total_time = 0.0
    for events in list(analog_events_by_qubit.values()) + list(virtual_events_by_qubit.values()):
        for event in events:
            total_time = max(total_time, event["start_us"] + event.get("width_us", 0.0))
    final_xlim_end = total_time * 1.05 if total_time > 0 else 1.0

    waveform_height = 1.8
    label_height = 0.5
    height_pattern = []
    for _ in qubit_ids:
        height_pattern.extend([label_height, waveform_height])
    fig_height = max(2.5, sum(height_pattern))
    fig = plt.figure(figsize=(14, fig_height))
    grid = fig.add_gridspec(nrows=len(height_pattern), ncols=1, height_ratios=height_pattern)

    colors = plt.cm.get_cmap("tab20")
    color_count = colors.N if hasattr(colors, "N") else 20

    last_wave_ax = None
    for idx, qubit in enumerate(qubit_ids):
        wave_ax = fig.add_subplot(grid[2 * idx + 1])
        label_ax = fig.add_subplot(grid[2 * idx])
        label_ax.axis("off")

        wave_ax.axhline(0.0, color="#cccccc", linewidth=0.8)
        wave_ax.set_ylabel(f"q{qubit}\nAmplitude")
        analog_events = analog_events_by_qubit.get(qubit, [])
        virt_events = virtual_events_by_qubit.get(qubit, [])
        rendered = []
        max_amp = 1e-3
        for event in analog_events:
            ts, i_samples, q_samples = sample_waveform(event, samples_per_us)
            times = event["start_us"] + ts
            peak = float(max(np.max(np.abs(i_samples)), np.max(np.abs(q_samples)) if np.any(q_samples) else 0.0))
            max_amp = max(max_amp, peak)
            rendered.append((event, times, i_samples, q_samples))

        for event_idx, (event, times, i_samples, q_samples) in enumerate(rendered):
            color = colors(event_idx % color_count)
            label_i = "I" if event_idx == 0 else ""
            label_q = "Q" if event_idx == 0 else ""
            wave_ax.plot(times, i_samples, color=color, linewidth=1.6, label=label_i)
            if np.any(q_samples):
                wave_ax.plot(times, q_samples, color=color, linewidth=1.2, linestyle="--", label=label_q)

        wave_ax.set_ylim(-max(0.2, 1.35 * max_amp), max(0.2, 1.2 * max_amp))
        wave_ax.set_xlim(0, final_xlim_end)

        for event in virt_events:
            x = event["start_us"]
            wave_ax.axvline(x, color="#8a2be2", linestyle="--", linewidth=1.0, alpha=0.7)

        labels = label_sequences.get(qubit, [])
        if labels:
            label_ax.set_xlim(0, 1)
            label_ax.set_ylim(0, 1)
            count = len(labels)
            if count == 1:
                x_positions = [0.5]
            else:
                x_positions = np.linspace(0.05, 0.95, count)
            for x_pos, entry in zip(x_positions, labels):
                if entry.get("virtual"):
                    text = f"rz({format_theta(entry.get('theta'))})"
                    color = "#4b0082"
                else:
                    if entry.get("gate") == "rx" and entry.get("theta") is not None:
                        text = f"rx({format_theta_compact(entry.get('theta'))})"
                    else:
                        text = entry.get("label", "pulse")
                    color = "#333333"
                label_ax.text(
                    x_pos,
                    0.5,
                    text,
                    rotation=90,
                    ha="center",
                    va="center",
                    fontsize=6,
                    color=color,
                    transform=label_ax.transAxes,
                )
            label_ax.set_xticks([])
            label_ax.set_yticks([])
        else:
            label_ax.set_xlim(0, 1)
            label_ax.set_ylim(0, 1)

        if idx == 0:
            handles, labels = wave_ax.get_legend_handles_labels()
            legend_labels = []
            legend_handles = []
            for handle, label in zip(handles, labels):
                if label and label not in legend_labels:
                    legend_labels.append(label)
                    legend_handles.append(handle)
            if legend_handles:
                wave_ax.legend(legend_handles, legend_labels, loc="upper right", frameon=False)

        wave_ax.tick_params(axis="x", which="both", labelbottom=False)
        last_wave_ax = wave_ax

    if last_wave_ax is not None:
        last_wave_ax.tick_params(axis="x", which="both", labelbottom=True)
        last_wave_ax.set_xlabel("Time (µs)")

    if title:
        fig.suptitle(title)
    fig.tight_layout(rect=(0, 0, 1, 0.98 if title else 1))
    fig.savefig(output_path, dpi=dpi)
    plt.close(fig)


def derive_output_path(pulses_path: Path, supplied_output: Optional[Path]):
    if supplied_output:
        return supplied_output
    return pulses_path.with_name(pulses_path.stem + "_plot.png")


def main(argv: Optional[Iterable[str]] = None):
    parser = argparse.ArgumentParser(description=__doc__ or "")
    parser.add_argument("circuit_pulses", type=Path, help="Circuit pulse JSON produced by QASMTrans")
    parser.add_argument("--output", type=Path, help="Output image path (defaults to <input>_plot.png)")
    parser.add_argument("--title", help="Optional plot title")
    parser.add_argument("--dpi", type=int, default=150, help="Figure DPI (default: 150)")
    parser.add_argument(
        "--samples-per-us",
        type=float,
        default=200.0,
        help="Sampling density for waveform rendering (default: 200 samples per microsecond)",
    )
    args = parser.parse_args(argv)

    pulse_doc = load_json(args.circuit_pulses)
    analog_events_by_qubit, virtual_events_by_qubit, label_sequences = build_qubit_events(pulse_doc)
    output_path = derive_output_path(args.circuit_pulses, args.output)
    plot_events(
        analog_events_by_qubit,
        virtual_events_by_qubit,
        label_sequences,
        args.title,
        output_path,
        args.dpi,
        args.samples_per_us,
    )
    print(f"Saved pulse timeline to {output_path}")


if __name__ == "__main__":
    main()
