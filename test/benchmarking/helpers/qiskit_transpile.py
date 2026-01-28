import os
import sys
import csv
import json
import time
import argparse
from qiskit import QuantumCircuit, transpile
from qiskit.transpiler import CouplingMap

def count_1q_2q_gates(circ):
    n1 = 0
    n2 = 0
    for inst, qargs, cargs in circ.data:
        if inst.name in ("barrier", "measure", "delay"):
            continue
        qn = len(qargs)
        if qn == 1:
            n1 += 1
        elif qn == 2:
            n2 += 1
    return n1, n2

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--qasm_dir", default="QiskitTest")
    parser.add_argument("--qasm_files", nargs="*", help="Explicit QASM files to process.")
    parser.add_argument("--csv_out", default="transpile_times.csv")
    parser.add_argument("--toronto_config", default="data/devices/ibmq_toronto.json")
    parser.add_argument("--brisbane_config", default="data/devices/ibm_brisbane.json")
    parser.add_argument("--toronto_threshold", type=int, default=27)
    parser.add_argument("--max_qubits", type=int, help="Skip circuits requiring more than this many qubits.")
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()

    row_fieldnames = [
        "circuit",
        "optimization_level",
        "time_ms",
        "depth",
        "one_qubit_gates",
        "two_qubit_gates",
        "backend",
        "backend_num_qubits",
        "circuit_qubits",
        "circuit_clbits",
    ]
    levels = [1, 2, 3]
    rows = []

    def load_device(path):
        if not os.path.isfile(path):
            raise FileNotFoundError(f"Device config not found: {path}")
        with open(path, "r") as fh:
            cfg = json.load(fh)

        coupling = cfg.get("cx_coupling") or cfg.get("coupling_map")
        if not coupling:
            raise ValueError(f"No coupling_map/cx_coupling found in device config {path}")

        edges = []
        for entry in coupling:
            if isinstance(entry, str):
                parts = entry.replace("-", "_").split("_")
            else:
                parts = entry
            if len(parts) != 2:
                raise ValueError(f"Unexpected coupling entry {entry!r} in {path}")
            u, v = map(int, parts)
            edges.append((u, v))

        cmap = CouplingMap(edges)
        basis = cfg.get("basis_gates") or ["rz", "sx", "x", "cx"]
        seen = set()
        basis_gates = []
        for gate in basis:
            if gate not in seen:
                basis_gates.append(gate)
                seen.add(gate)

        return {
            "map": cmap,
            "basis": basis_gates,
            "name": cfg.get("name") or os.path.splitext(os.path.basename(path))[0],
            "num_qubits": int(cfg.get("num_qubits", cmap.size())),
        }

    devices = {
        "toronto": load_device(args.toronto_config),
        "brisbane": load_device(args.brisbane_config),
    }

    qasm_paths = []
    if args.qasm_files:
        for path in args.qasm_files:
            if not os.path.isfile(path):
                raise FileNotFoundError(f"QASM file not found: {path}")
            qasm_paths.append(os.path.abspath(path))
    else:
        if not os.path.isdir(args.qasm_dir):
            raise FileNotFoundError(f"Directory not found: {args.qasm_dir}")
        for fname in sorted(os.listdir(args.qasm_dir)):
            if fname.endswith(".qasm"):
                qasm_paths.append(os.path.abspath(os.path.join(args.qasm_dir, fname)))

    for idx, path in enumerate(qasm_paths, start=1):
        fname = os.path.basename(path)
        device_info = None
        try:
            qc = QuantumCircuit.from_qasm_file(path)
            n_qubits = qc.num_qubits
            if args.max_qubits is not None and n_qubits > args.max_qubits:
                continue
            n_clbits = qc.num_clbits
            device_key = "toronto" if n_qubits <= args.toronto_threshold else "brisbane"
            device_info = devices[device_key]
            print(f"[Qiskit] ({idx}/{len(qasm_paths)}) {fname}: {n_qubits} qubits, using {device_info['name']}")
            sys.stdout.flush()
            for lvl in levels:
                t0 = time.perf_counter()
                tc = transpile(
                    qc,
                    coupling_map=device_info["map"],
                    basis_gates=device_info["basis"],
                    optimization_level=lvl,
                    seed_transpiler=args.seed
                )
                dt_ms = (time.perf_counter() - t0) * 1000.0
                depth = tc.depth()
                n1, n2 = count_1q_2q_gates(tc)
                print(f"  level {lvl}: {dt_ms:.3f} ms, depth={depth}, 1q={n1}, 2q={n2}")
                sys.stdout.flush()
                rows.append({
                    "circuit": fname,
                    "optimization_level": lvl,
                    "time_ms": f"{dt_ms:.3f}",
                    "depth": depth,
                    "one_qubit_gates": n1,
                    "two_qubit_gates": n2,
                    "backend": device_info["name"],
                    "backend_num_qubits": device_info["num_qubits"],
                    "circuit_qubits": n_qubits,
                    "circuit_clbits": n_clbits
                })
        except Exception as e:
            rows.append({
                "circuit": fname,
                "optimization_level": "error",
                "time_ms": str(e),
                "depth": "",
                "one_qubit_gates": "",
                "two_qubit_gates": "",
                "backend": device_info["name"] if device_info else "",
                "backend_num_qubits": device_info["num_qubits"] if device_info else "",
                "circuit_qubits": "",
                "circuit_clbits": ""
            })

    with open(args.csv_out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=row_fieldnames)
        writer.writeheader()
        writer.writerows(rows)

if __name__ == "__main__":
    main()
