#!/usr/bin/env bash
# =============================================================================
# 004_core_biom_Pycnandra_V2.sh
# Decontamination par PYC-NEG, retrait de tous les ASV observes dans PYC-NEG,
# exclusion du controle de la table biologique, rarefaction, core metrics et
# exports TSV/BIOM. Les echantillons -10- ne sont pas presents depuis 001.
# =============================================================================
set -Eeuo pipefail
shopt -s nullglob
IFS=$'\n\t'
PROJECT_NAME="Pycnandra_V2"
PROJECT_DIR="/nvme/bio/data_fungi/${PROJECT_NAME}"
RESULTS_DIR="${PROJECT_DIR}/02_amplicon_pipeline"
TMPDIR_BASE="${PROJECT_DIR}/tmp"
QIIME_ENV="qiime2-amplicon-2025.7"
BIOM_ENV="biom-format"
MARKERS=("16S" "ITS")
# Mettre a une profondeur choisie apres lecture de 002. 0 = profondeur auto :
# minimum de reads parmi les 12 echantillons biologiques decontamines.
SAMPLING_DEPTH_16S=0
SAMPLING_DEPTH_ITS=0
RUN_DECONTAM=true
RUN_RAREFY=true
RUN_CORE_METRICS=true
RUN_EXPORT=true

eval "$(conda shell.bash hook)"
mkdir -p "$TMPDIR_BASE"; export TMPDIR="$TMPDIR_BASE"
qiime_run(){ conda run -n "$QIIME_ENV" qiime "$@"; }
biom_run(){ conda run -n "$BIOM_ENV" biom "$@"; }
require(){ [[ -f "$1" ]] || { echo "Artefact absent: $1" >&2; exit 1; }; }
for marker in "${MARKERS[@]}"; do
  lc="$(tr '[:upper:]' '[:lower:]' <<< "$marker")"
  MDIR="${RESULTS_DIR}/${lc}"; CORE="${MDIR}/05_qiime2/core"; TREE="${MDIR}/05_qiime2/tree/rooted-tree.qza"; META="${MDIR}/04_database_files/sample-metadata_${lc}.tsv"; DECON="${MDIR}/05_qiime2/decontam"; EXP="${MDIR}/05_qiime2/export"; LOG="${RESULTS_DIR}/logs/004_core_biom_${lc}.log"
  mkdir -p "$DECON" "$EXP"
  INPUT="$CORE/table.qza"; require "$INPUT"; require "$TREE"; require "$META"
  NEG_IDS="$DECON/negative_asv_ids.qza"; NO_NEG_ASV="$DECON/table_no_negative_asvs.qza"; BIO="$DECON/table_no_negative_asvs_biological_only.qza"; REMOVED="$DECON/asvs_removed_from_PYC_NEG.tsv"
  if [[ "$RUN_DECONTAM" == true ]]; then
  rm -f "$NO_NEG_ASV" "$BIO" \
        "$DECON/table_PYC_NEG_only.qza" \
        "$DECON/negative_asv_ids.txt"
  rm -rf "$DECON/export_neg"

  # 1. Isoler l'échantillon PYC-NEG dans la table DADA2.
  qiime_run feature-table filter-samples \
    --i-table "$INPUT" \
    --m-metadata-file "$META" \
    --p-where "[sample_type]='negative_control'" \
    --o-filtered-table "$DECON/table_PYC_NEG_only.qza"

  # 2. Exporter sa table et récupérer les identifiants des ASV présents.
  qiime_run tools export \
    --input-path "$DECON/table_PYC_NEG_only.qza" \
    --output-path "$DECON/export_neg"

  biom_run convert \
    -i "$DECON/export_neg/feature-table.biom" \
    -o "$DECON/negative_control_table.tsv" \
    --to-tsv

  awk 'NR > 2 && $1 != "" {print $1}' \
    "$DECON/negative_control_table.tsv" \
    > "$DECON/negative_asv_ids.txt"

  [[ -s "$DECON/negative_asv_ids.txt" ]] || {
    echo "Aucun ASV détecté dans PYC-NEG pour $marker : arrêt pour vérification." >&2
    exit 1
  }

  # 3. Exclure ces ASV de la table entière.
  qiime_run feature-table filter-features \
    --i-table "$INPUT" \
    --m-metadata-file "$DECON/negative_asv_ids.txt" \
    --p-exclude-ids \
    --o-filtered-table "$NO_NEG_ASV"

  # 4. Ne conserver que les échantillons biologiques pour la diversité.
  qiime_run feature-table filter-samples \
    --i-table "$NO_NEG_ASV" \
    --m-metadata-file "$META" \
    --p-where "[sample_type]='biological'" \
    --o-filtered-table "$BIO"

  cp "$DECON/negative_asv_ids.txt" "$REMOVED"
fi
  require "$BIO"
  summary_qzv="$DECON/table_biological_decontam_summary.qzv"; summary_dir="$DECON/table_biological_decontam_summary"
  rm -f "$summary_qzv"; rm -rf "$summary_dir"; qiime_run feature-table summarize --i-table "$BIO" --m-sample-metadata-file "$META" --o-visualization "$summary_qzv"; qiime_run tools export --input-path "$summary_qzv" --output-path "$summary_dir"
  auto_depth=$(python - "$summary_dir/data/sample-frequency-detail.csv" <<'PY'
import pandas as pd,sys
x=pd.read_csv(sys.argv[1]); print(int(x.iloc[:,1].min()))
PY
)
  case "$marker" in 16S) configured="$SAMPLING_DEPTH_16S";; ITS) configured="$SAMPLING_DEPTH_ITS";; esac
  depth="$auto_depth"; [[ "$configured" -gt 0 ]] && depth="$configured"
  RAR="$DECON/RarTable_depth${depth}.qza"; CM="$DECON/core-metrics-depth${depth}"; EXPCM="$EXP/core-metrics-depth${depth}"
  if [[ "$RUN_RAREFY" == true ]]; then rm -f "$RAR"; qiime_run feature-table rarefy --i-table "$BIO" --p-sampling-depth "$depth" --p-no-with-replacement --o-rarefied-table "$RAR"; fi
  require "$RAR"
  if [[ "$RUN_CORE_METRICS" == true ]]; then rm -rf "$CM"; qiime_run diversity core-metrics-phylogenetic --i-phylogeny "$TREE" --i-table "$BIO" --p-sampling-depth "$depth" --m-metadata-file "$META" --p-n-jobs-or-threads 1 --output-dir "$CM"; fi
  if [[ "$RUN_EXPORT" == true ]]; then
    rm -rf "$EXPCM"; mkdir -p "$EXPCM/rarefied_table"
    qiime_run tools export --input-path "$RAR" --output-path "$EXPCM/rarefied_table"
    biom_run convert -i "$EXPCM/rarefied_table/feature-table.biom" -o "$EXPCM/rarefied_table/table-from-biom.tsv" --to-tsv
    sed '1d; s/^#OTU ID/ASV_ID/' "$EXPCM/rarefied_table/table-from-biom.tsv" > "$EXPCM/rarefied_table/ASV.tsv"
    for a in faith_pd_vector shannon_vector evenness_vector observed_features_vector bray_curtis_distance_matrix jaccard_distance_matrix unweighted_unifrac_distance_matrix weighted_unifrac_distance_matrix bray_curtis_pcoa_results jaccard_pcoa_results unweighted_unifrac_pcoa_results weighted_unifrac_pcoa_results; do qiime_run tools export --input-path "$CM/${a}.qza" --output-path "$EXPCM/$a"; done
    for a in bray_curtis_emperor jaccard_emperor unweighted_unifrac_emperor weighted_unifrac_emperor; do qiime_run tools export --input-path "$CM/${a}.qzv" --output-path "$EXPCM/$a"; done
  fi
  echo "$marker : decontamination PYC-NEG, 12 echantillons biologiques, profondeur=${depth}, resultats=${EXPCM}"
done
