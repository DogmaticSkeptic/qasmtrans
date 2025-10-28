#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
BUILD_DIR="$REPO_ROOT/build"

DEVICE_NAME="ibm_brisbane"
DEVICE_JSON="../data/devices/ibm_brisbane.json"
DEVICE_QUBITS=127
COUNTS=(1 2 4 6)

NWQSIM_ROOT="$REPO_ROOT/../NWQ-Sim"
NWQSIM_BUILD="$NWQSIM_ROOT/build"
DEFAULT_NWQSIM_EXE="$NWQSIM_BUILD/qasm/nwq_qasm"
NWQSIM_EXE="${NWQSIM_EXE:-$DEFAULT_NWQSIM_EXE}"
NWQSIM_SHOTS="${NWQSIM_SHOTS:-1024}"
NWQSIM_BACKEND="${NWQSIM_BACKEND:-CPU}"
MAPOMATIC_LIMIT="${MAPOMATIC_LIMIT:-50000}"
declare -a NWQSIM_EXTRA_ARGS_ARRAY=()
if [ -n "${NWQSIM_EXTRA_ARGS:-}" ]; then
    NWQSIM_EXTRA_ARGS_STR="${NWQSIM_EXTRA_ARGS}"
    IFS=$' \t\n' read -r -a NWQSIM_EXTRA_ARGS_ARRAY <<<"$NWQSIM_EXTRA_ARGS_STR"
    unset NWQSIM_EXTRA_ARGS_STR
fi

ADAPTIVE_DIR="$REPO_ROOT/../AdaptiveData"
OUTPUT_CSV="$REPO_ROOT/adaptive_timing_ibm_brisbane.csv"
TMP_DIR="$REPO_ROOT/tmp"

if [ ! -d "$ADAPTIVE_DIR" ]; then
    echo "AdaptiveData directory not found at $ADAPTIVE_DIR" >&2
    exit 1
fi

if [ ! -x "$NWQSIM_EXE" ]; then
    echo "NWQ-Sim executable not found at $NWQSIM_EXE. Build NWQ-Sim or set NWQSIM_EXE." >&2
    exit 1
fi

mkdir -p "$(dirname "$OUTPUT_CSV")"
mkdir -p "$TMP_DIR"
echo "Circuit,Qubit Number,Device,Device Qubit number,Number of circuits compiled to device,Time to transpile single circuit to device (ms),Time to transpile circuit number to device (ms),Circuit partitioning time (ms),Average Fidelity,Minimum Fidelity,Maximum Fidelity" >"$OUTPUT_CSV"

cd "$BUILD_DIR"

python_qubit_count() {
    python3 - "$1" <<'PY'
import re, sys
path = sys.argv[1]
count = None
with open(path, 'r') as fh:
    for line in fh:
        m = re.search(r'qreg\s+\w+\[(\d+)\]', line)
        if m:
            count = int(m.group(1))
            break
if count is None:
    raise SystemExit(f"Could not determine qubit count in {path}")
print(count)
PY
}

now_seconds() {
    python3 - <<'PY'
import time
print("{:.9f}".format(time.perf_counter()))
PY
}

elapsed_ms() {
    python3 - "$1" "$2" <<'PY'
import sys
start = float(sys.argv[1])
end = float(sys.argv[2])
print(int((end - start) * 1000))
PY
}

abspath() {
    python3 - "$1" <<'PY'
import os
import sys
print(os.path.abspath(sys.argv[1]))
PY
}

compute_fidelity_stats() {
    local subchip_dir="$1"
    local count="$2"
    if [ ! -d "$subchip_dir" ]; then
        echo ",,"
        return
    fi

    shopt -s nullglob
    local -a fidelities=()
    local subchip_json
    for subchip_json in "$subchip_dir"/circuit*_subchip.json; do
        local base="${subchip_json%.json}"
        local circuit_qasm="${base}.qasm"
        if [ ! -f "$circuit_qasm" ]; then
            continue
        fi
        local -a nwq_cmd=("$NWQSIM_EXE" --backend "$NWQSIM_BACKEND" --shots "$NWQSIM_SHOTS" --sim dm \
            --device "$subchip_json" -q "$circuit_qasm" --fidelity)
        if [ "${#NWQSIM_EXTRA_ARGS_ARRAY[@]}" -gt 0 ]; then
            nwq_cmd+=("${NWQSIM_EXTRA_ARGS_ARRAY[@]}")
        fi
        local nwq_output
        if ! nwq_output="$("${nwq_cmd[@]}" 2>&1)"; then
            printf '%s\n' "$nwq_output" >&2
            echo "NWQ-Sim simulation failed for $circuit_qasm" >&2
            exit 1
        fi
        local fidelity
        if ! fidelity="$(printf '%s\n' "$nwq_output" | python3 -c '
import re
import sys
text = sys.stdin.read()
match = re.search(r"State Fidelity:\s*([0-9.eE+-]+)", text)
if not match:
    sys.exit(1)
print(match.group(1))
')"; then
            printf '%s\n' "$nwq_output" >&2
            echo "Failed to extract fidelity from NWQ-Sim output for $circuit_qasm" >&2
            exit 1
        fi
        fidelities+=("$fidelity")
    done
    shopt -u nullglob

    if [ "${#fidelities[@]}" -eq 0 ]; then
        echo ",,"
        return
    fi

    python3 - "$count" "${fidelities[@]}" <<'PY'
import sys
count = int(sys.argv[1])
values = [float(v) for v in sys.argv[2:]]
avg = sum(values) / len(values)
if count > 1:
    lo = min(values)
    hi = max(values)
    print(f"{avg:.6f},{lo:.6f},{hi:.6f}")
else:
    print(f"{avg:.6f},,")
PY
}

for circuit_path in "$ADAPTIVE_DIR"/*.qasm; do
    [ -e "$circuit_path" ] || continue
    circuit_name="$(basename "$circuit_path")"
    circuit_base="${circuit_name%.*}"
    qubit_count="$(python_qubit_count "$circuit_path")"

    single_ms=""

    for count in "${COUNTS[@]}"; do
        inputs=()
        for ((i = 0; i < count; ++i)); do
            inputs+=("-i" "$circuit_path")
        done

        timestamp="$(date +"%Y%m%d_%H%M%S")"
        output_base="../data/output_qasm_file/${circuit_base}_N${count}_${timestamp}"
        log_file="$(mktemp "$TMP_DIR/qasmtrans.XXXXXX")"

        start_time="$(now_seconds)"
        if ! ./qasmtrans "${inputs[@]}" -m ibmq -c "$DEVICE_JSON" -o "${output_base}.qasm" -full_fidelity -mapomatic_limit "$MAPOMATIC_LIMIT" >"$log_file" 2>&1; then
            cat "$log_file" >&2
            rm -f "$log_file"
            echo "Transpilation failed for $circuit_name (count=$count)" >&2
            exit 1
        fi
        end_time="$(now_seconds)"

        runtime_ms="$(elapsed_ms "$start_time" "$end_time")"
        partition_ms="$(grep -o 'partition_ms=[0-9]*' "$log_file" | tail -1 | cut -d= -f2 || echo "")"
        total_ms="$(grep -o 'total_ms=[0-9]*' "$log_file" | tail -1 | cut -d= -f2 || echo "")"

        if [ -z "$partition_ms" ]; then
            partition_ms=0
        fi
        if [ -z "$total_ms" ]; then
            total_ms="$runtime_ms"
        fi

        if [ -z "$single_ms" ]; then
            single_ms="$total_ms"
        fi

        output_abs="$(abspath "$output_base")"
        fidelity_fields="$(compute_fidelity_stats "${output_abs}_subchips" "$count")"
        IFS=',' read -r avg_fidelity min_fidelity max_fidelity <<<"$fidelity_fields"

        printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
            "$circuit_name" \
            "$qubit_count" \
            "$DEVICE_NAME" \
            "$DEVICE_QUBITS" \
            "$count" \
            "$single_ms" \
            "$total_ms" \
            "$partition_ms" \
            "$avg_fidelity" \
            "$min_fidelity" \
            "$max_fidelity" >>"$OUTPUT_CSV"

        rm -f "$log_file"
    done
done

echo "CSV written to $OUTPUT_CSV"
