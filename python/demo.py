#!/usr/bin/env python3
"""Showcase the qasmtrans_core Python API (transpilation + pulse plotting)."""
from __future__ import annotations

import json
from pathlib import Path
from tempfile import NamedTemporaryFile

REPO = Path(__file__).resolve().parents[1]

import qasmtrans as qt  # type: ignore
import qasmtrans.plot_pulses as plot_pulses  # type: ignore


def demo_transpile_from_file():
    qasm_path = REPO / "data" / "test_benchmark" / "adder_n4.qasm"
    backend_path = REPO / "data" / "devices" / "ibmq_toronto.json"

    opts = qt.TranspileOptions()
    opts.backend_config = str(backend_path)
    opts.output_path = str(REPO / "data" / "output_demo")
    opts.verbose = 1
    opts.disable_mapomatic = False
    opts.mapomatic_limit = 200

    result = qt.transpile_qasm(str(qasm_path), opts)
    print("\n=== Transpile from file ===")
    print("Log:", result.log)
    print("QASM path:", result.output_qasm_path)
    print("First few QASM lines:")
    for line in result.output_qasm.splitlines()[:8]:
        print("  ", line)
    return result


def demo_transpile_from_string():
    qasm_path = REPO / "data" / "test_benchmark" / "adder_n4.qasm"
    backend_path = REPO / "data" / "devices" / "ibmq_toronto.json"
    qasm_text = qasm_path.read_text()

    opts = qt.TranspileOptions()
    opts.backend_config = str(backend_path)
    opts.output_path = str(REPO / "data" / "output_demo_string")
    opts.disable_mapomatic = True  # show how to toggle Mapomatic
    opts.verbose = 0

    result = qt.transpile_qasm(qasm_text, opts)
    print("\n=== Transpile from raw string ===")
    print("Log:", result.log)
    print("QASM length:", len(result.output_qasm))
    return result


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


def build_minimal_pulse_template(path: Path):
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


def demo_pulse_plot():
    # Small 2-qubit circuit to generate pulses from the minimal template.
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
    opts.output_path = str(REPO / "data" / "output_demo_pulses")
    opts.limited_qubits = True

    result = qt.transpile_qasm(qasm_text, opts)
    print("\n=== Pulse compilation + plot ===")
    print("Log:", result.log)
    print("Pulse JSON path:", result.pulse_path)

    pulse_doc = result.pulse_doc if result.pulse_doc is not None else json.loads(result.pulse_schedule or "{}")
    analog_events, virtual_events, labels = plot_pulses.build_qubit_events(pulse_doc)
    plot_path = REPO / "data" / "output_demo_pulses.png"
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


def main():
    demo_transpile_from_file()
    demo_transpile_from_string()
    demo_pulse_plot()
    print("\nNote: QICK hardware demo is in python/demo_qick.py")


if __name__ == "__main__":
    main()
