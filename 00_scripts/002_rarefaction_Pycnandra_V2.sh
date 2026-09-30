#!/usr/bin/env bash
# =============================================================================
# 002_rarefaction_Pycnandra_V2.sh
# Courbes de rarefaction alpha apres retrait des ASV presents dans PYC-NEG.
# Les profondeurs maximales sont determinees automatiquement a partir de la
# table decontaminee, et non reprises du projet Araucaria.
# =============================================================================
set -Eeuo pipefail
PROJECT_NAME="Pycnandra_V2"
PROJECT_DIR="/nvme/bio/data_fungi/${PROJECT_NAME}"
RESULTS_DIR="${PROJECT_DIR}/02_amplicon_pipeline"
QIIME_ENV="qiime2-amplicon-2025.7"
MARKERS=("16S" "ITS")
eval "$(conda shell.bash hook)"
conda activate "$QIIME_ENV"
for marker in "${MARKERS[@]}"; do
  lc="$(tr '[:upper:]' '[:lower:]' <<< "$marker")"
  QIIME_DIR="${RESULTS_DIR}/${lc}/05_qiime2"
  DATABASE="${RESULTS_DIR}/${lc}/04_database_files"
  TABLE="${QIIME_DIR}/decontam/table_no_negative_asvs_biological_only.qza"
  TREE="${QIIME_DIR}/tree/rooted-tree.qza"
  META="${DATABASE}/sample-metadata_${lc}.tsv"
  OUT="${QIIME_DIR}/visual/alpha-rarefaction_${lc}_decontam_biological_only.qzv"
  [[ -f "$TABLE" && -f "$TREE" && -f "$META" ]] || { echo "Artefact manquant pour $marker. Lancez d'abord 004." >&2; exit 1; }
  max_depth=$(qiime feature-table summarize --i-table "$TABLE" --o-visualization "${QIIME_DIR}/visual/_tmp_table_summary_${lc}.qzv" >/dev/null 2>&1; qiime tools export --input-path "${QIIME_DIR}/visual/_tmp_table_summary_${lc}.qzv" --output-path "${QIIME_DIR}/visual/_tmp_table_summary_${lc}" >/dev/null 2>&1; python - "${QIIME_DIR}/visual/_tmp_table_summary_${lc}/data/sample-frequency-detail.csv" <<'PY'
import pandas as pd,sys
x=pd.read_csv(sys.argv[1]); print(int(x.iloc[:,1].max()))
PY
)
  rm -rf "${QIIME_DIR}/visual/_tmp_table_summary_${lc}"
  qiime diversity alpha-rarefaction --i-table "$TABLE" --i-phylogeny "$TREE" --p-max-depth "$max_depth" --p-min-depth 1 --m-metadata-file "$META" --o-visualization "$OUT"
  echo "$marker : rarefaction alpha ecrite dans $OUT (max-depth=$max_depth)"
done
