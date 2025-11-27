#!/usr/bin/env bash
set -euo pipefail

# usage: ./meta_lrc_nine.sh path/to/lrc_for_nine_runners.cpp k p1 p2 ... [output.txt]
# Example: ./meta_lrc_nine.sh lrc_for_nine_runners.cpp 8 53 59 results_nine.txt

if [ "$#" -lt 3 ]; then
  echo "Usage: $0 path/to/lrc_for_nine_runners.cpp k p1 p2 ... [output.txt]"
  exit 1
fi

SRC="$1"
K="$2"
shift 2

# Optional output file: use last argument if it ends in .txt
OUTFILE="lrc_nine_runs.txt"
if [ "$#" -ge 1 ] && [[ "${!#}" == *.txt ]]; then
  OUTFILE="${!#}"
  set -- "${@:1:$(($#-1))}"  # remove last argument
fi

if [ "$#" -lt 1 ]; then
  echo "No primes provided."
  exit 1
fi

PRIMES=("$@")
WORKBIN="./lrc_nine_tmp"
COMPILE_ERR="compile_error.tmp"
RUN_LOG="run_output.tmp"

# Create or truncate the output file
: > "$OUTFILE"

for p in "${PRIMES[@]}"; do
  ts=$(date -Iseconds)
  echo "=== RUN START: $ts | prime=$p | k=$K ===" >> "$OUTFILE"
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
  start_time=$(date +%s)
  start_iso=$(date -Iseconds)
  echo "Execution start: $start_iso" >> "$OUTFILE"

  if "$WORKBIN" >> "$OUTFILE" 2>&1; then
    status="OK"
  else
    status="NONZERO_EXIT"
  fi

  end_time=$(date +%s)
  end_iso=$(date -Iseconds)
  runtime=$((end_time - start_time))

  echo "Execution end: $end_iso" >> "$OUTFILE"
  echo "Status: $status" >> "$OUTFILE"
  echo "Runtime_seconds: ${runtime}" >> "$OUTFILE"
  echo "=== RUN END: prime=$p | runtime=${runtime}s ===" >> "$OUTFILE"
  echo >> "$OUTFILE"
done

# Cleanup temp files
rm -f "$WORKBIN" "$COMPILE_ERR" "$RUN_LOG" || true

echo "All runs complete. Combined output written to $OUTFILE"