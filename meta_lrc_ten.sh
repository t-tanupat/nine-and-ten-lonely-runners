#!/usr/bin/env bash
set -euo pipefail

# usage: ./meta_lrc_ten.sh path/to/lrc_for_ten_runners.cpp K p1 p2 p3 ... [output.txt]
# Example: ./meta_lrc_ten.sh lrc_for_ten_runners.cpp 9 137 139 149 results_ten.txt

if [ "$#" -lt 3 ]; then
  echo "Usage: $0 path/to/lrc_for_ten_runners.cpp K p1 p2 ... [output.txt]"
  exit 1
fi

SRC="$1"
K="$2"
shift 2

# Optional output file
OUTFILE="lrc_ten_runs.txt"
if [ "$#" -ge 1 ] && [[ "${!#}" == *.txt ]]; then
  OUTFILE="${!#}"
  set -- "${@:1:$(($#-1))}"
fi

if [ "$#" -lt 1 ]; then
  echo "No primes provided."
  exit 1
fi

PRIMES=("$@")
WORKBIN="./lrc_ten_tmp"
COMPILE_ERR="compile_error.tmp"

# Create or truncate output file
: > "$OUTFILE"

for p in "${PRIMES[@]}"; do
  ts=$(date -Iseconds)
  echo "=== RUN START: $ts | prime=$p | K=$K ===" >> "$OUTFILE"
  echo "Compiling (PRIME=$p K=$K) ..." >> "$OUTFILE"

  if ! g++ -std=c++17 -O2 -DPRIME="$p" -DK="$K" "$SRC" -o "$WORKBIN" 2> "$COMPILE_ERR"; then
    echo "=== COMPILE FAILED for prime=$p ===" >> "$OUTFILE"
    echo "Compiler output:" >> "$OUTFILE"
    sed -n '1,200p' "$COMPILE_ERR" >> "$OUTFILE"
    echo "=== END RUN: $ts | prime=$p | compile failed ===" >> "$OUTFILE"
    echo >> "$OUTFILE"
    continue
  fi

  echo "=== EXECUTING: binary for PRIME=$p ===" >> "$OUTFILE"
  start_time_ns=$(date +%s%N)
  start_iso=$(date -Iseconds)
  echo "Execution start: $start_iso" >> "$OUTFILE"

  if "$WORKBIN" >> "$OUTFILE" 2>&1; then
    status="OK"
  else
    status="NONZERO_EXIT"
  fi

  end_time_ns=$(date +%s%N)
  end_iso=$(date -Iseconds)
  elapsed_ns=$((end_time_ns - start_time_ns))
  elapsed_s=$(awk -v n="$elapsed_ns" 'BEGIN{ printf("%.3f", n/1e9) }')

  echo "Execution end: $end_iso" >> "$OUTFILE"
  echo "Status: $status" >> "$OUTFILE"
  echo "Runtime_seconds: ${elapsed_s}" >> "$OUTFILE"
  echo "=== RUN END: prime=$p | runtime=${elapsed_s}s ===" >> "$OUTFILE"
  echo >> "$OUTFILE"
done

# Cleanup
rm -f "$WORKBIN" "$COMPILE_ERR" || true

echo "All runs complete. Combined output written to $OUTFILE"