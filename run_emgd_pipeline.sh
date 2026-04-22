#!/usr/bin/env bash
set -euo pipefail

DEPENDENCY_CSV="data/processed/gene_dependency_train_example.csv"
PROTEOMICS_CSV="data/processed/proteomics_train.csv"

GROUP_DIR="results/emgd_groups"
MARKER_DIR="results/emgd_markers"
DVALUE_DIR="results/d_values"

mkdir -p "$GROUP_DIR" "$MARKER_DIR" "$DVALUE_DIR"

echo "Step 1/3: Building EMGD groups"
py scripts/02_build_emgd_groups.py \
  --dependency "$DEPENDENCY_CSV" \
  --proteomics "$PROTEOMICS_CSV" \
  --output-dir "$GROUP_DIR"

for strategy in median quartile high_low extremes cluster; do
  echo "Step 2/3: Generating markers using ${strategy} strategy"
  py scripts/03_generate_emgd_protein_markers.py \
    --dependency "$DEPENDENCY_CSV" \
    --proteomics "$PROTEOMICS_CSV" \
    --groups "$GROUP_DIR/emgd_${strategy}_groups.json" \
    --output "$MARKER_DIR/emgd_${strategy}_marker_db.csv"

  echo "Step 3/3: Calculating D-values for ${strategy}"
  py scripts/04_calculate_d_values.py \
    --markers "$MARKER_DIR/emgd_${strategy}_marker_db.csv" \
    --proteomics "$PROTEOMICS_CSV" \
    --output "$DVALUE_DIR/emgd_${strategy}_dvalues.csv" \
    --missing-at-random
done

echo "Pipeline completed successfully."