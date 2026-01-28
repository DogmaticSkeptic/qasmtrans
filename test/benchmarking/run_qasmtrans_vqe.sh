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
OUTPUT_CSV="${OUTPUT_CSV:-$REPO_ROOT/data/vqe_uccsd_qasmtrans_1q_opt_compare.csv}"
OPTIMIZE_2Q_CANCEL="${OPTIMIZE_2Q_CANCEL:-0}"
OPTIMIZE_COMMUTE_2Q="${OPTIMIZE_COMMUTE_2Q:-0}"
OPTIMIZE_2Q_SYNTH="${OPTIMIZE_2Q_SYNTH:-0}"

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

echo "circuit,base_reported_ms,base_elapsed_ms,base_one_qubit,base_two_qubit,base_depth,opt_reported_ms,opt_elapsed_ms,opt_one_qubit,opt_two_qubit,opt_depth" >"$OUTPUT_CSV"

run_qasmtrans() {
    local label="$1"
    local output_qasm="$2"
    shift 2
    local extra_flags=("$@")

    local cmd=(
        "$QASMTRANS_BIN"
        -i "$circuit_path"
        -m ibmq
        -c "$DEVICE_JSON"
        -o "$output_qasm"
        -v 1
        --disable_mapomatic
    )

    if [ "$label" = "opt" ]; then
        if [ "$OPTIMIZE_2Q_CANCEL" = "1" ]; then
            cmd+=(--optimize-2q-cancel)
        fi
        if [ "$OPTIMIZE_COMMUTE_2Q" = "1" ]; then
            cmd+=(--optimize-commute-2q)
        fi
        if [ "$OPTIMIZE_2Q_SYNTH" = "1" ]; then
            cmd+=(--optimize-2q-synth)
        fi
    fi

    if [ "${#extra_flags[@]}" -gt 0 ]; then
        cmd+=("${extra_flags[@]}")
    fi

    echo "  $label: ${cmd[*]}" >&2
    local start_ns
    start_ns="$($PYTHON_BIN - <<'PY'
import time
print(int(time.perf_counter() * 1e9))
PY
)"

    local output
    if ! output="$("${cmd[@]}" 2>&1 | tee "$TMP_DIR/${circuit_name}_${label}_attempt${attempt}.log")"; then
        echo "    QASMTrans exited with error" >&2
        return 1
    fi

    local end_ns
    end_ns="$($PYTHON_BIN - <<'PY'
import time
print(int(time.perf_counter() * 1e9))
PY
)"

    local elapsed_ms
    elapsed_ms="$($PYTHON_BIN - <<PY
start_ns=int("${start_ns}")
end_ns=int("${end_ns}")
print((end_ns - start_ns)/1e6)
PY
)"

    local reported_line
    reported_line="$(echo "$output" | grep '^ total QASMTrans time' | head -n 1)"
    local reported_ms
    if [ -n "$reported_line" ]; then
        reported_ms="$($PYTHON_BIN - "$reported_line" <<'PY'
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

    if [ ! -f "$output_qasm" ]; then
        echo "    Output QASM not found at $output_qasm" >&2
        return 1
    fi
    local metrics_values
    if ! metrics_values="$($PYTHON_BIN - "$output_qasm" <<'PY'
import sys
from qiskit import QuantumCircuit

path = sys.argv[1]
qc = QuantumCircuit.from_qasm_file(path)
one = two = 0
for inst, qargs, _ in qc.data:
    if inst.name in {"barrier", "measure", "delay"}:
        continue
    if len(qargs) == 1:
        one += 1
    elif len(qargs) == 2:
        two += 1
print(f"{one}:{two}:{qc.depth()}")
PY
)"; then
        echo "    Failed to parse QASM output for metrics" >&2
        return 1
    fi

    local one_qubit two_qubit depth
    IFS=':' read -r one_qubit two_qubit depth <<<"$metrics_values"

    local formatted_reported formatted_elapsed
    formatted_reported="$($PYTHON_BIN - "$reported_ms" <<'PY'
import sys
val = float(sys.argv[1])
print(f"{val:.3f}")
PY
)"
    formatted_elapsed="$($PYTHON_BIN - "$elapsed_ms" <<'PY'
import sys
val = float(sys.argv[1])
print(f"{val:.3f}")
PY
)"

    echo "${formatted_reported}:${formatted_elapsed}:${one_qubit}:${two_qubit}:${depth}"
}

for circuit_path in "$INPUT_DIR"/*.qasm; do
    [ -e "$circuit_path" ] || continue
    circuit_name="$(basename "$circuit_path")"
    echo "=== $circuit_name ==="

    success=0
    for ((attempt = 1; attempt <= ATTEMPTS; ++attempt)); do
        base_output_qasm="$TMP_DIR/${circuit_name%.qasm}_base_attempt${attempt}.qasm"
        opt_output_qasm="$TMP_DIR/${circuit_name%.qasm}_opt_attempt${attempt}.qasm"

        base_result="$(run_qasmtrans "base" "$base_output_qasm")" || continue
        opt_result="$(run_qasmtrans "opt" "$opt_output_qasm" --optimize-1q)" || continue

        IFS=':' read -r base_reported base_elapsed base_one base_two base_depth <<<"$base_result"
        IFS=':' read -r opt_reported opt_elapsed opt_one opt_two opt_depth <<<"$opt_result"

        echo "    Success: base=${base_reported} ms, opt=${opt_reported} ms"
        echo "${circuit_name},${base_reported},${base_elapsed},${base_one},${base_two},${base_depth},${opt_reported},${opt_elapsed},${opt_one},${opt_two},${opt_depth}" >>"$OUTPUT_CSV"
        success=1
        break
    done

    if [ "$success" -eq 0 ]; then
        echo "    Failed after ${ATTEMPTS} attempt(s); leaving blank entry"
        echo "${circuit_name},,,,,,,,,," >>"$OUTPUT_CSV"
    fi
done

echo "Saved metrics to ${OUTPUT_CSV}"

"$PYTHON_BIN" - "$OUTPUT_CSV" <<'PY'
import csv
import statistics
import sys

path = sys.argv[1]
with open(path, newline="") as fh:
    rows = list(csv.DictReader(fh))

def collect(key):
    vals = []
    for row in rows:
        val = row.get(key, "")
        if val == "":
            continue
        vals.append(float(val))
    return vals

def report(label, values):
    if not values:
        print(f"{label}: n=0")
        return
    mean = statistics.mean(values)
    print(f"{label}: n={len(values)} mean={mean:.3f}")

report("QASMTrans base elapsed ms", collect("base_elapsed_ms"))
report("QASMTrans opt elapsed ms", collect("opt_elapsed_ms"))
report("QASMTrans base depth", collect("base_depth"))
report("QASMTrans opt depth", collect("opt_depth"))
PY
