#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python}"

INPUT_DIR="${INPUT_DIR:-$REPO_ROOT/data/test_benchmark/vqe_uccsd_generated}"
DEVICE_JSON="${DEVICE_JSON:-$REPO_ROOT/data/devices/ibm_brisbane.json}"
QASMTRANS_BIN="${QASMTRANS_BIN:-$REPO_ROOT/build/QASMTrans}"
TIMEOUT_MS="${TIMEOUT_MS:-1000}"
ATTEMPTS="${ATTEMPTS:-3}"
OUTPUT_CSV="${OUTPUT_CSV:-$REPO_ROOT/data/vqe_uccsd_qasmtrans_metrics.csv}"

if [ ! -x "$QASMTRANS_BIN" ]; then
    echo "QASMTrans binary not found at $QASMTRANS_BIN" >&2
    exit 1
fi

if [ ! -d "$INPUT_DIR" ]; then
    echo "Input directory not found: $INPUT_DIR" >&2
    exit 1
fi

TMP_DIR="$(mktemp -d "$REPO_ROOT/tmp/qasmtrans_vqe.XXXXXX")"
trap 'rm -rf "$TMP_DIR"' EXIT

echo "circuit,qasmtrans_reported_ms,qasmtrans_elapsed_ms,qiskit_o1_time_ms,qiskit_o2_time_ms,qiskit_o3_time_ms,one_qubit,two_qubit,depth" >"$OUTPUT_CSV"

for circuit_path in "$INPUT_DIR"/*.qasm; do
    [ -e "$circuit_path" ] || continue
    circuit_name="$(basename "$circuit_path")"
    echo "=== $circuit_name ==="

    success=0
    for ((attempt = 1; attempt <= ATTEMPTS; ++attempt)); do
        output_qasm="$TMP_DIR/${circuit_name%.qasm}_attempt${attempt}.qasm"
        cmd=(
            "$QASMTRANS_BIN"
            -i "$circuit_path"
            -m ibmq
            -c "$DEVICE_JSON"
            -o "$output_qasm"
            -v 1
            --disable_mapomatic
        )
        echo "  Attempt $attempt: ${cmd[*]}"
        start_ns="$("$PYTHON_BIN" - <<'PY'
import time
print(int(time.perf_counter() * 1e9))
PY
)"

        if ! output="$("${cmd[@]}" 2>&1 | tee "$TMP_DIR/${circuit_name}_attempt${attempt}.log")"; then
            echo "    QASMTrans exited with error"
            continue
        fi

        end_ns="$("$PYTHON_BIN" - <<'PY'
import time
print(int(time.perf_counter() * 1e9))
PY
)"
        elapsed_ms="$("$PYTHON_BIN" - <<PY
start_ns=int("${start_ns}")
end_ns=int("${end_ns}")
print((end_ns - start_ns)/1e6)
PY
)"

        metrics_line="$(echo "$output" | grep '^\[metrics\]')"
        if [ -z "$metrics_line" ]; then
            echo "    Metrics not found, retrying"
            continue
        fi

        reported_line="$(echo "$output" | grep '^ total QASMTrans time')"
        if [ -n "$reported_line" ]; then
            reported_ms="$("$PYTHON_BIN" - "$reported_line" <<'PY'
import re
import sys
line = sys.argv[1]
match = re.search(r"([-+]?[0-9]*\.?[0-9]+)\s*ms", line)
print(match.group(1) if match else "")
PY
)"
            if [ -z "$reported_ms" ]; then
                reported_ms="$elapsed_ms"
            fi
        else
            reported_ms="$elapsed_ms"
        fi

        metrics_values="$("$PYTHON_BIN" - "$metrics_line" <<'PY'
import re
import sys
line = sys.argv[1]
keys = ("one_qubit_gates", "two_qubit_gates", "depth")
values = []
for key in keys:
    match = re.search(rf"{key}=(\d+)", line)
    values.append(match.group(1) if match else "")
print(":".join(values))
PY
)"
        OLD_IFS="$IFS"
        IFS=':' read -r one_qubit two_qubit depth <<<"$metrics_values"
        IFS="$OLD_IFS"

        formatted_reported="$("$PYTHON_BIN" - "$reported_ms" <<'PY'
import sys
val = float(sys.argv[1])
print(f"{val:.3f}")
PY
)"
        formatted_elapsed="$("$PYTHON_BIN" - "$elapsed_ms" <<'PY'
import sys
val = float(sys.argv[1])
print(f"{val:.3f}")
PY
)"

        qiskit_times="$("$PYTHON_BIN" - "$circuit_path" "$DEVICE_JSON" <<'PY'
import json
import sys
import time
from qiskit import QuantumCircuit, transpile
from qiskit.transpiler import CouplingMap

circuit_path, device_json = sys.argv[1:3]
with open(device_json, "r") as fh:
    cfg = json.load(fh)
coupling = cfg.get("cx_coupling") or cfg.get("coupling_map")
edges = []
for entry in coupling:
    if isinstance(entry, str):
        parts = entry.replace("-", "_").split("_")
        u, v = map(int, parts)
    else:
        u, v = entry
    edges.append((u, v))
cmap = CouplingMap(edges)
basis = cfg.get("basis_gates") or ["rz", "sx", "x", "cx"]

qc = QuantumCircuit.from_qasm_file(circuit_path)
times = []
for level in (1, 2, 3):
    start = time.perf_counter()
    transpile(qc, coupling_map=cmap, basis_gates=basis, optimization_level=level)
    elapsed = (time.perf_counter() - start) * 1e3
    times.append(f"{elapsed:.3f}")
print(":".join(times))
PY
)"
        OLD_IFS="$IFS"
        IFS=':' read -r qiskit_o1 qiskit_o2 qiskit_o3 <<<"$qiskit_times"
        IFS="$OLD_IFS"

        echo "    Success: reported=${formatted_reported} ms, elapsed=${formatted_elapsed} ms, qiskit(o1/o2/o3)=(${qiskit_o1}, ${qiskit_o2}, ${qiskit_o3}) ms"
        echo "${circuit_name},${formatted_reported},${formatted_elapsed},${qiskit_o1},${qiskit_o2},${qiskit_o3},${one_qubit},${two_qubit},${depth}" >>"$OUTPUT_CSV"
        success=1
        break
    done

    if [ "$success" -eq 0 ]; then
        echo "    Failed after ${ATTEMPTS} attempt(s); leaving blank entry"
        echo "${circuit_name},,,,,,,," >>"$OUTPUT_CSV"
    fi
done

echo "Saved metrics to ${OUTPUT_CSV}"
