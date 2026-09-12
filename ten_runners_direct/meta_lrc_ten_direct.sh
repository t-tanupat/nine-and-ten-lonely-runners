#!/usr/bin/env bash
set -euo pipefail

# Usage: ./meta_lrc_ten_direct.sh path/to/lrc_for_ten_runners_direct.cpp K p1 p2 ... [output.txt]
# Example: ./meta_lrc_ten_direct.sh lrc_for_ten_runners_direct.cpp 9 137 139 149 results_ten_direct.txt

if [ "$#" -lt 3 ]; then
  echo "Usage: $0 path/to/lrc_for_ten_runners_direct.cpp K p1 p2 ... [output.txt]"
  exit 1
fi

SRC="$1"
K="$2"
shift 2

OUTFILE="lrc_ten_direct_runs.txt"
if [[ "${!#}" == *.txt ]]; then
  OUTFILE="${!#}"
  set -- "${@:1:$(($#-1))}"
fi

if [ "$#" -lt 1 ]; then
  echo "No primes provided."
  exit 1
fi

PRIMES=("$@")
WORKDIR=$(mktemp -d "${TMPDIR:-/tmp}/lrc_ten_direct.XXXXXX")
WORKBIN="$WORKDIR/lrc_ten_direct_tmp"
COMPILE_ERR="$WORKDIR/compile_error.tmp"
trap 'rm -f "$WORKBIN" "$COMPILE_ERR"; rmdir "$WORKDIR"' EXIT

# Create or truncate output file.
: > "$OUTFILE"

for p in "${PRIMES[@]}"; do
  ts=$(date '+%Y-%m-%dT%H:%M:%S%z')
  echo "=== RUN START: $ts | prime=$p | K=$K ===" >> "$OUTFILE"
  echo "Compiling (PRIME=$p K=$K) ..." >> "$OUTFILE"

  if ! g++ -std=c++17 -O2 -pthread -DPRIME="$p" -DK="$K" "$SRC" -o "$WORKBIN" 2> "$COMPILE_ERR"; then
    echo "=== COMPILE FAILED for prime=$p ===" >> "$OUTFILE"
    echo "Compiler output:" >> "$OUTFILE"
    sed -n '1,200p' "$COMPILE_ERR" >> "$OUTFILE"
    echo "=== END RUN: $ts | prime=$p | compile failed ===" >> "$OUTFILE"
    echo >> "$OUTFILE"
    continue
  fi

  echo "=== EXECUTING: binary for PRIME=$p ===" >> "$OUTFILE"
  # BSD date on macOS does not support nanoseconds (%N).
  start_time=$(date +%s)
  start_iso=$(date '+%Y-%m-%dT%H:%M:%S%z')
  echo "Execution start: $start_iso" >> "$OUTFILE"

  if "$WORKBIN" >> "$OUTFILE" 2>&1; then
    status="OK"
  else
    status="NONZERO_EXIT"
  fi

  end_time=$(date +%s)
  end_iso=$(date '+%Y-%m-%dT%H:%M:%S%z')
  elapsed_s=$((end_time - start_time))

  echo "Execution end: $end_iso" >> "$OUTFILE"
  echo "Status: $status" >> "$OUTFILE"
  echo "Runtime_seconds: ${elapsed_s}" >> "$OUTFILE"
  echo "=== RUN END: prime=$p | runtime=${elapsed_s}s ===" >> "$OUTFILE"
  echo >> "$OUTFILE"
done

echo "All runs complete. Combined output written to $OUTFILE"
