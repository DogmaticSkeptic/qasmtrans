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
    events_by_qubit: Dict[str, List[dict]] = defaultdict(list)
    for item in pulse_doc.get("schedule", []):
        pulse_id = item.get("pulse_id")
        pulse = library.get(pulse_id, {})
        qubits = item.get("qubits") or []
        if not qubits:
            continue
        start = float(item.get("start_time", 0.0))
        duration = float(item.get("duration", pulse.get("width", 0.0)))
        if duration <= 0:
            continue
        amplitude = float(pulse.get("amplitude", 0.0))
        label = item.get("gate") or pulse.get("gate") or pulse_id or "pulse"
        shape = (pulse.get("shape") or "const").lower()
        for qubit in qubits:
            events_by_qubit[str(qubit)].append(
                {
                    "start_us": start * 1e6,
                    "width_us": duration * 1e6,
                    "amplitude": amplitude,
                    "shape": shape,
                    "waveform_type": (pulse.get("waveform_type") or shape),
                    "label": f"{label}\n({pulse_id})" if pulse_id else label,
                    "samples_i": [float(v) for v in pulse.get("samples_i", [])],
                    "samples_q": [float(v) for v in pulse.get("samples_q", [])],
                }
            )
    for events in events_by_qubit.values():
        events.sort(key=lambda entry: entry["start_us"])
    return dict(events_by_qubit)


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


def plot_events(events_by_qubit, title: Optional[str], output_path: Path, dpi: int, samples_per_us: float):
    if not events_by_qubit:
        raise ValueError("No pulse events were found in the schedule.")

    qubit_ids = sorted(events_by_qubit.keys(), key=lambda q: int(q) if str(q).isdigit() else str(q))
    total_time = max(
        (event["start_us"] + event["width_us"] for events in events_by_qubit.values() for event in events),
        default=0.0,
    )

    fig_height = max(2.5, 2.0 * len(qubit_ids))
    fig, axes = plt.subplots(len(qubit_ids), 1, sharex=True, figsize=(14, fig_height))
    if len(qubit_ids) == 1:
        axes = [axes]

    colors = plt.cm.get_cmap("tab20")
    color_count = colors.N if hasattr(colors, "N") else 20

    for idx, qubit in enumerate(qubit_ids):
        ax = axes[idx]
        ax.axhline(0.0, color="#cccccc", linewidth=0.8)
        ax.set_ylabel(f"q{qubit}\nAmplitude")
        events = events_by_qubit[qubit]
        rendered = []
        max_amp = 1e-3
        for event in events:
            ts, i_samples, q_samples = sample_waveform(event, samples_per_us)
            times = event["start_us"] + ts
            peak = float(max(np.max(np.abs(i_samples)), np.max(np.abs(q_samples)) if np.any(q_samples) else 0.0))
            max_amp = max(max_amp, peak)
            rendered.append((event, times, i_samples, q_samples))

        for event_idx, (event, times, i_samples, q_samples) in enumerate(rendered):
            color = colors(event_idx % color_count)
            label_i = "I" if event_idx == 0 else ""
            label_q = "Q" if event_idx == 0 else ""
            ax.plot(times, i_samples, color=color, linewidth=1.6, label=label_i)
            if np.any(q_samples):
                ax.plot(times, q_samples, color=color, linewidth=1.2, linestyle="--", label=label_q)
            ax.text(
                event["start_us"] + event["width_us"] / 2.0,
                1.05 * max_amp,
                event["label"],
                ha="center",
                va="bottom",
                fontsize=8,
            )
        ax.set_ylim(-1.3 * max_amp, 1.3 * max_amp)
        if idx == 0:
            handles, labels = ax.get_legend_handles_labels()
            legend_labels = []
            legend_handles = []
            for handle, label in zip(handles, labels):
                if label and label not in legend_labels:
                    legend_labels.append(label)
                    legend_handles.append(handle)
            if legend_handles:
                ax.legend(legend_handles, legend_labels, loc="upper right", frameon=False)

    axes[-1].set_xlabel("Time (µs)")
    axes[-1].set_xlim(0, total_time * 1.05 if total_time > 0 else 1.0)
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
    events_by_qubit = build_qubit_events(pulse_doc)
    output_path = derive_output_path(args.circuit_pulses, args.output)
    plot_events(events_by_qubit, args.title, output_path, args.dpi, args.samples_per_us)
    print(f"Saved pulse timeline to {output_path}")


if __name__ == "__main__":
    main()
