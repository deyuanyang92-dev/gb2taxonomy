#!/usr/bin/env bash
# Reproduce the benchmark: download the taxa from NCBI, run g2t with and without voucher
# reconciliation, time each run and score it against the published tables.
#
#   bash benchmark/run_benchmark.sh WORKDIR [NCBI_EMAIL]
#
# Needs network access (NCBI E-utilities). Set NCBI_API_KEY for faster downloads.
# Set SUPERCRUNCH_PARSE_LOCI=/path/to/SuperCRUNCH/supercrunch-scripts/Parse_Loci.py to add SuperCRUNCH's parser.
# Results: WORKDIR/results/<dataset>.json, WORKDIR/results/timing.tsv (appended to when it exists)
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
work="$(realpath -m "${1:?usage: run_benchmark.sh WORKDIR [NCBI_EMAIL]}")"
[ -n "${2:-}" ] && export NCBI_EMAIL="$2"
mkdir -p "$work/results"
gold="$here/gold/gold_standard.tsv"
timing="$work/results/timing.tsv"
[ -s "$timing" ] || printf "dataset\tmode\tseconds\tmax_rss_mb\n" > "$timing"

timed() {  # timed NAME MODE CMD...
  local name=$1 mode=$2; shift 2
  local t0=$(date +%s.%N)
  /usr/bin/time -f "%M" -o "$work/results/.rss" "$@"
  local t1=$(date +%s.%N)
  printf "%s\t%s\t%.1f\t%.0f\n" "$name" "$mode" "$(echo "$t1 - $t0" | bc)" "$(echo "$(cat "$work/results/.rss") / 1024" | bc)" >> "$timing"
}

# dataset name | gold dataset | input
#   Haplotaxidae: the 347 records cited by Martin et al. (2024), bundled in examples/
#   others: every marker record of the taxon in NCBI nuccore (g2t-download default selection)
declare -A INPUT=(
  [haplotaxidae]="$here/../examples/haplotaxidae/haplotaxidae.gb"
  [hesionidae]="Hesionidae"
  [opheliidae]="Opheliidae"
  [phyllodocida]="Phyllodocida"
  [priapulidae]="Priapulidae"
  [nymphaeaceae]="Nymphaeaceae"
  [ginkgoaceae]="Ginkgoaceae"
  [ophioglossaceae]="Ophioglossaceae"
)
# plant families have no published specimen table: they show gene-type coverage, specimens and run time
declare -A GOLD=([haplotaxidae]=Haplotaxidae [hesionidae]=Hesionidae [opheliidae]=Opheliidae [phyllodocida]=Polynoidae
  [priapulidae]=none [nymphaeaceae]=none [ginkgoaceae]=none [ophioglossaceae]=none)

# DATASETS="priapulidae nymphaeaceae" re-runs only those (timing rows are appended)
for name in ${DATASETS:-haplotaxidae hesionidae opheliidae phyllodocida priapulidae nymphaeaceae ginkgoaceae ophioglossaceae}; do
  src="${INPUT[$name]}"
  if [ -f "$src" ]; then
    in="$src"
  else
    timed "$name" download g2t-download -t "$src" -o "$work/gb_$name" -q
    in="$work/gb_$name/markers"
  fi
  timed "$name" g2t g2t -i "$in" -o "$work/out_$name" --stream --quiet
  timed "$name" g2t_exact g2t -i "$in" -o "$work/out_${name}_exact" --stream --quiet --skip_reconcile
  timed "$name" g2t_high g2t -i "$in" -o "$work/out_${name}_high" --stream --quiet --reconcile_min_confidence high
  python "$here/evaluate.py" --gold "$gold" --dataset "${GOLD[$name]}" --bootstrap 1000 \
    --run "$work/out_$name" --run-exact "$work/out_${name}_exact" --run-high "$work/out_${name}_high" \
    --details "$work/results/${name}_genes.tsv" ${SUPERCRUNCH_PARSE_LOCI:+--supercrunch "$SUPERCRUNCH_PARSE_LOCI"} \
    > "$work/results/$name.json"
done
rm -f "$work/results/.rss"
python "$here/summarize.py" "$work"
