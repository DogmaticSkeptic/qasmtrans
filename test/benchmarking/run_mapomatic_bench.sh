#!/usr/bin/env bash
set -euo pipefail

CIRCUITS=(vqe8 qpe9 sat7 shor7)
MODES=(nomap product full hybrid)
INPUT_DIR="../tmp"
DEVICE_JSON="../data/devices/ibm_brisbane.json"
OUTPUT_DIR="../data/output_qasm_file"
CSV_OUT="../data/mapomatic_benchmarks.csv"
NWQ_BIN="../../NWQ-Sim/build/qasm/nwq_qasm"
MAPOMATIC_LIMIT="${MAPOMATIC_LIMIT:-50000}"

# helper to timestamp (ms)
timestamp_ms() {
  python3 - <<'PY'
import time
print(int(time.time()*1000))
PY
}

# ensure output directory exists
mkdir -p "$OUTPUT_DIR"

echo "circuit,mode,total_ms,mapomatic_ms,fidelity" > "$CSV_OUT.tmp"

for circuit in "${CIRCUITS[@]}"; do
  input_qasm="$INPUT_DIR/${circuit}.qasm"
  if [[ ! -f "$input_qasm" ]]; then
    echo "Missing input $input_qasm" >&2
    exit 1
  fi
  for mode in "${MODES[@]}"; do
    out_prefix="$OUTPUT_DIR/${circuit}_${mode}"
    mapomatic_ms=0
    case "$mode" in
      nomap)
        cmd=(./qasmtrans -i "$input_qasm" -m ibmq -c "$DEVICE_JSON" -o "${out_prefix}.qasm" --disable_mapomatic)
        ;;
      product)
        cmd=(./qasmtrans -i "$input_qasm" -m ibmq -c "$DEVICE_JSON" -o "${out_prefix}.qasm" -mapomatic_limit "$MAPOMATIC_LIMIT")
        ;;
      full)
        cmd=(./qasmtrans -i "$input_qasm" -m ibmq -c "$DEVICE_JSON" -o "${out_prefix}.qasm" -full_fidelity -mapomatic_limit "$MAPOMATIC_LIMIT")
        ;;
      hybrid)
        cmd=(./qasmtrans -i "$input_qasm" -m ibmq -c "$DEVICE_JSON" -o "${out_prefix}.qasm" -cp_mode hybrid -mapomatic_limit "$MAPOMATIC_LIMIT")
        ;;
      *)
        echo "Unknown mode $mode" >&2
        exit 1
        ;;
    esac

    if [[ "$mode" != "nomap" ]]; then
      cmd+=(-v 2)
    fi

    start_ms=$(timestamp_ms)
    # capture output for timing stats
    map_output=$("${cmd[@]}" 2>&1)
    end_ms=$(timestamp_ms)
    total_ms=$((end_ms - start_ms))

    # parse mapomatic timing if present
    if [[ "$mode" != "nomap" ]]; then
      if echo "$map_output" | grep -qi "Mapomatic timing"; then
        timing_line=$(echo "$map_output" | grep -i "Mapomatic timing" | head -n1)
        prep=$(echo "$timing_line" | awk -F'prep ' '{print $2}' | awk -F',' '{print $1}')
        embedding=$(echo "$timing_line" | awk -F'embedding ' '{print $2}' | awk -F',' '{print $1}')
        scoring=$(echo "$timing_line" | awk -F'scoring ' '{print $2}' | awk -F',' '{print $1}')
        apply=$(echo "$timing_line" | awk -F'apply ' '{print $2}' | awk -F',' '{print $1}')
        mapomatic_ms=$(python3 - <<PY
prep=float("${prep:-0}")
embed=float("${embedding:-0}")
scor=float("${scoring:-0}")
apply=float("${apply:-0}")
print(f"{prep+embed+scor+apply:.3f}")
PY
)
      fi
    fi

    # run NWQ-Sim
    subchip_dir="${out_prefix}_subchips"
    device_json="$subchip_dir/circuit00_subchip.json"
    qasm_file="$subchip_dir/circuit00_subchip.qasm"
    if [[ ! -f "$device_json" ]]; then
      echo "Missing device JSON $device_json" >&2
      exit 1
    fi
    fidelity_line=$($NWQ_BIN --backend CPU --shots 1024 --sim dm --device "$device_json" --fidelity -q "$qasm_file" | tail -n 1)
    fidelity=$(echo "$fidelity_line" | awk '{print $NF}')

    echo "$circuit,$mode,$total_ms,$mapomatic_ms,$fidelity" >> "$CSV_OUT.tmp"
  done

done

mv "$CSV_OUT.tmp" "$CSV_OUT"
echo "Results written to $CSV_OUT"
