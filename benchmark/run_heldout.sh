#!/usr/bin/env bash
# Held-out evaluation: four published tables not used while the rules were developed (rules frozen at tag
# frozen-heldout-v1). Each table is scored inside a realistic background download of its clade.
#
#   bash benchmark/run_heldout.sh WORKDIR
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
work="$(realpath -m "${1:?usage: run_heldout.sh WORKDIR}")"
mkdir -p "$work/results"
gold="$here/gold/heldout_gold.tsv"
timing="$work/results/timing.tsv"
printf "dataset\tmode\tseconds\tmax_rss_mb\n" > "$timing"
timed() {
  local name=$1 mode=$2; shift 2
  local t0=$(date +%s.%N)
  /usr/bin/time -f "%M" -o "$work/results/.rss" "$@"
  local t1=$(date +%s.%N)
  printf "%s\t%s\t%.1f\t%.0f\n" "$name" "$mode" "$(echo "$t1 - $t0" | bc)" "$(echo "$(cat "$work/results/.rss") / 1024" | bc)" >> "$timing"
}
# name | taxon | extra query
run() {
  local name=$1 taxon=$2 extra=${3:-}
  [ "${OFFLINE:-0}" = 1 ] || timed "$name" download g2t-download -t "$taxon" -o "$work/gb_$name" -q ${extra:+--query "$extra"} --tag bg
  local in="$work/gb_$name/bg"
  timed "$name" g2t g2t -i "$in" -o "$work/out_$name" --stream --quiet
  timed "$name" g2t_exact g2t -i "$in" -o "$work/out_${name}_exact" --stream --quiet --skip_reconcile
  timed "$name" g2t_high g2t -i "$in" -o "$work/out_${name}_high" --stream --quiet --reconcile_min_confidence high
}
run fionidae Cladobranchia
run nipponacmea Lottiidae "all[filter] NOT txid225164[Organism]"
run leptochilus Polypodiaceae
run angraecum Angraecinae
rm -f "$work/results/.rss"
