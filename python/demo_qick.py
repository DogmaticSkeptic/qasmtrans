#!/usr/bin/env python3
"""Demonstrate pulse generation, plotting, and QICK emission via qasmtrans_core."""
from __future__ import annotations

import json
from pathlib import Path
from tempfile import NamedTemporaryFile

REPO = Path(__file__).resolve().parents[1]

import qasmtrans as qt  # type: ignore
import qasmtrans.plot_pulses as plot_pulses  # type: ignore


def build_minimal_pulse_template(path: Path):
    """Write a minimal pulse template covering x, rz, and cx on qubits 0/1."""
    template = {
        "name": "demo_template",
        "version": "0.1",
        "pulse_definitions": [
            {"id": "rz_q0", "gate": "rz", "qubits": [0], "shape": "virtual", "waveform_type": "virtual", "width": 0.0, "amplitude": 0.0, "virtual": True},
            {"id": "rz_q1", "gate": "rz", "qubits": [1], "shape": "virtual", "waveform_type": "virtual", "width": 0.0, "amplitude": 0.0, "virtual": True},
            {"id": "x_q0", "gate": "x", "qubits": [0], "shape": "const", "waveform_type": "const", "width": 1e-6, "amplitude": 1.0},
            {"id": "x_q1", "gate": "x", "qubits": [1], "shape": "const", "waveform_type": "const", "width": 1e-6, "amplitude": 1.0},
            {"id": "cx_q0_q1", "gate": "cx", "qubits": [0, 1], "shape": "const", "waveform_type": "const", "width": 2e-6, "amplitude": 1.0},
        ],
    }
    path.write_text(json.dumps(template, indent=2))
    return path


def build_minimal_backend(path: Path):
    backend = {
        "name": "demo_backend",
        "version": "0.1",
        "num_qubits": 2,
        "basis_gates": ["x", "rz", "cx"],
        "gate_lens": {"x0": 1.0, "x1": 1.0, "rz0": 0.0, "rz1": 0.0, "cx0_1": 2.0},
        "gate_errs": {"x0": 0.0, "x1": 0.0, "rz0": 0.0, "rz1": 0.0, "cx0_1": 0.0},
        "cx_coupling": ["0_1"],
    }
    path.write_text(json.dumps(backend, indent=2))
    return path


def demo_qick():
    # Simple 2-qubit QASM to match the template definitions
    qasm_text = """OPENQASM 2.0;
include "qelib1.inc";
qreg q[2];
creg c[2];
x q[0];
cx q[0],q[1];
rz(1.5708) q[1];
measure q[0] -> c[0];
measure q[1] -> c[1];
"""

    with NamedTemporaryFile("w", suffix=".json", delete=False) as tmp_backend, \
            NamedTemporaryFile("w", suffix=".json", delete=False) as tmp_template:
        backend_path = build_minimal_backend(Path(tmp_backend.name))
        pulse_template_path = build_minimal_pulse_template(Path(tmp_template.name))

    opts = qt.TranspileOptions()
    opts.backend_config = str(backend_path)
    opts.pulse_template_path = str(pulse_template_path)
    opts.output_path = str(REPO / "data" / "output_demo_qick")
    opts.allow_parameterized_merge = True
    opts.limited_qubits = True

    result = qt.transpile_qasm(qasm_text, opts)
    print("Log:", result.log)
    print("Pulse path:", result.pulse_path)

    # Plot pulses to a PNG (no interactive window).
    pulse_doc = result.pulse_doc if result.pulse_doc is not None else json.loads(result.pulse_schedule or "{}")
    analog_events, virtual_events, labels = plot_pulses.build_qubit_events(pulse_doc)
    plot_path = REPO / "data" / "output_demo_qick_pulses.png"
    plot_pulses.plot_events(
        analog_events_by_qubit=analog_events,
        virtual_events_by_qubit=virtual_events,
        label_sequences=labels,
        title="QASMTrans Pulse Schedule",
        output_path=plot_path,
        dpi=150,
        samples_per_us=200.0,
    )
    print(f"Pulse plot saved to {plot_path}")

    pulse_doc = result.pulse_doc if result.pulse_doc is not None else json.loads(result.pulse_schedule or "{}")
    qick_cfg = {"qubit_gen_map": {"0": 0, "1": 1}, "pulse_freq": 6.0e9}
    try:
        summary = qt.emit_qick(pulse_doc, qick_cfg, run_enabled=False, summary_only=True)
        print("emit_qick summary (no hardware run):", json.dumps(summary, indent=2))
    except ImportError as exc:
        print("QICK not available (install qick on a QICK board):", exc)
    except Exception as exc:
        print("emit_qick failed:", exc)


if __name__ == "__main__":
    demo_qick()
