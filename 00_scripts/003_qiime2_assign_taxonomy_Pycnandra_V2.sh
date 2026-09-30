#!/usr/bin/env bash
# =============================================================================
# 003_qiime2_assign_taxonomy_Pycnandra_V2.sh
# Taxonomie des ASV Pycnandra_V2 : SILVA 16S V4 et UNITE ITS7/ITS4.
# La taxonomie est assignee a tous les ASV, y compris ceux des controles,
# afin que 004 puisse auditer les ASV retires via PYC-NEG.
# =============================================================================
set -Eeuo pipefail
shopt -s nullglob
IFS=$'\n\t'
trap 'rc=$?; echo "[ERREUR] ligne ${LINENO}: ${BASH_COMMAND}" >&2; exit "$rc"' ERR
PROJECT_NAME="Pycnandra_V2"
PROJECT_DIR="/nvme/bio/data_fungi/${PROJECT_NAME}"
RESULTS_DIR="${PROJECT_DIR}/02_amplicon_pipeline"
TMPDIR_BASE="${PROJECT_DIR}/tmp"
QIIME2_ENV="qiime2-amplicon-2025.7"
THREADS=8
SILVA_VERSION="138.2"
SILVA_TARGET="SSURef_NR99"
UNITE_VERSION="10.0_2025-02-19"
UNITE_RELEASE_DATE="2025-02-19"
UNITE_TAXON_GROUP="fungi"
UNITE_CLUSTER_ID="dynamic"
UNITE_SINGLETONS=true
PRIMER_F_16S="GTGCCAGCMGCCGCGGTAA"
PRIMER_R_16S="GGACTACHVGGGTWTCTAAT"
PRIMER_F_ITS="GTGARTCATCGAATCTTTG"
PRIMER_R_ITS="TCCTCCGCTTATTGATATGC"
PRIMER_IDENTITY=0.80
MIN_LENGTH_16S=100
MIN_LENGTH_ITS=50
TAXONOMY_CONFIDENCE=0.70
READS_PER_BATCH="auto"
REBUILD_16S_CLASSIFIER=true
REBUILD_ITS_CLASSIFIER=true
DOWNLOAD_16S_SILVA=true
DOWNLOAD_ITS_UNITE=true

eval "$(conda shell.bash hook)"
conda activate "$QIIME2_ENV"
mkdir -p "$TMPDIR_BASE"; export TMPDIR="$TMPDIR_BASE"
require(){ [[ -s "$1" ]] || { echo "Fichier absent/vide: $1" >&2; exit 1; }; }
setup(){
  marker="$1"; lc="$(tr '[:upper:]' '[:lower:]' <<< "$marker")"
  MARKER_DIR="${RESULTS_DIR}/${lc}"; CORE="${MARKER_DIR}/05_qiime2/core"; META="${MARKER_DIR}/04_database_files/sample-metadata_${lc}.tsv"
  TAX_DIR="${MARKER_DIR}/05_qiime2/taxonomy"; DB="${TAX_DIR}/database"; CLS="${TAX_DIR}/classifier"; RES="${TAX_DIR}/results"; EXP="${MARKER_DIR}/05_qiime2/export/taxonomy"
  TABLE="${CORE}/table.qza"; REP="${CORE}/rep-seqs.qza"; mkdir -p "$DB" "$CLS" "$RES" "$EXP"; require "$TABLE"; require "$REP"; require "$META"
}
outputs(){
  prefix="$1"; tax="$2"
  qiime metadata tabulate --m-input-file "$tax" --o-visualization "${RES}/${prefix}_taxonomy.qzv"
  qiime taxa barplot --i-table "$TABLE" --i-taxonomy "$tax" --m-metadata-file "$META" --o-visualization "${RES}/${prefix}_taxa-barplot.qzv"
  rm -rf "${EXP}/${prefix}_taxonomy"; qiime tools export --input-path "$tax" --output-path "${EXP}/${prefix}_taxonomy"
}
# 16S
setup 16S
SILVA_RNA="${DB}/silva-${SILVA_VERSION}-ssu-rna.qza"; SILVA_DNA="${DB}/silva-${SILVA_VERSION}-ssu-dna.qza"; SILVA_TAX="${DB}/silva-${SILVA_VERSION}-tax.qza"; SILVA_EXTRACT="${DB}/silva-${SILVA_VERSION}-515F-806R.qza"; SILVA_UNIQ="${DB}/silva-${SILVA_VERSION}-515F-806R-uniq.qza"; SILVA_UTAX="${DB}/silva-${SILVA_VERSION}-515F-806R-uniq-tax.qza"; SILVA_CLS="${CLS}/silva-${SILVA_VERSION}-515F-806R-classifier.qza"; SILVA_OUT="${RES}/16S_SILVA_${SILVA_VERSION}_taxonomy.qza"
if [[ "$DOWNLOAD_16S_SILVA" == true ]]; then qiime rescript get-silva-data --p-version "$SILVA_VERSION" --p-target "$SILVA_TARGET" --o-silva-sequences "$SILVA_RNA" --o-silva-taxonomy "$SILVA_TAX"; qiime rescript reverse-transcribe --i-rna-sequences "$SILVA_RNA" --o-dna-sequences "$SILVA_DNA"; fi
if [[ "$REBUILD_16S_CLASSIFIER" == true ]]; then qiime feature-classifier extract-reads --i-sequences "$SILVA_DNA" --p-f-primer "$PRIMER_F_16S" --p-r-primer "$PRIMER_R_16S" --p-identity "$PRIMER_IDENTITY" --p-min-length "$MIN_LENGTH_16S" --p-n-jobs "$THREADS" --p-read-orientation forward --o-reads "$SILVA_EXTRACT"; qiime rescript dereplicate --i-sequences "$SILVA_EXTRACT" --i-taxa "$SILVA_TAX" --p-mode uniq --o-dereplicated-sequences "$SILVA_UNIQ" --o-dereplicated-taxa "$SILVA_UTAX"; qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads "$SILVA_UNIQ" --i-reference-taxonomy "$SILVA_UTAX" --o-classifier "$SILVA_CLS"; fi
qiime feature-classifier classify-sklearn --i-reads "$REP" --i-classifier "$SILVA_CLS" --p-n-jobs "$THREADS" --p-reads-per-batch "$READS_PER_BATCH" --p-confidence "$TAXONOMY_CONFIDENCE" --p-read-orientation auto --o-classification "$SILVA_OUT"; outputs "16S_SILVA_${SILVA_VERSION}" "$SILVA_OUT"
# ITS
setup ITS
USEQ="${DB}/unite-${UNITE_VERSION}-seqs.qza"; UTAX="${DB}/unite-${UNITE_VERSION}-tax.qza"; UEX="${DB}/unite-${UNITE_VERSION}-ITS7-ITS4.qza"; UUNIQ="${DB}/unite-${UNITE_VERSION}-ITS7-ITS4-uniq.qza"; UU_TAX="${DB}/unite-${UNITE_VERSION}-ITS7-ITS4-uniq-tax.qza"; UCLS="${CLS}/unite-${UNITE_VERSION}-ITS7-ITS4-classifier.qza"; UOUT="${RES}/ITS_UNITE_${UNITE_VERSION}_taxonomy.qza"
if [[ "$DOWNLOAD_ITS_UNITE" == true ]]; then opt=(--p-no-singletons); [[ "$UNITE_SINGLETONS" == true ]] && opt=(--p-singletons); qiime rescript get-unite-data --p-version "$UNITE_RELEASE_DATE" --p-taxon-group "$UNITE_TAXON_GROUP" --p-cluster-id "$UNITE_CLUSTER_ID" "${opt[@]}" --o-sequences "$USEQ" --o-taxonomy "$UTAX"; fi
if [[ "$REBUILD_ITS_CLASSIFIER" == true ]]; then qiime feature-classifier extract-reads --i-sequences "$USEQ" --p-f-primer "$PRIMER_F_ITS" --p-r-primer "$PRIMER_R_ITS" --p-identity "$PRIMER_IDENTITY" --p-min-length "$MIN_LENGTH_ITS" --p-n-jobs "$THREADS" --p-read-orientation both --o-reads "$UEX"; qiime rescript dereplicate --i-sequences "$UEX" --i-taxa "$UTAX" --p-mode uniq --o-dereplicated-sequences "$UUNIQ" --o-dereplicated-taxa "$UU_TAX"; qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads "$UUNIQ" --i-reference-taxonomy "$UU_TAX" --o-classifier "$UCLS"; fi
qiime feature-classifier classify-sklearn --i-reads "$REP" --i-classifier "$UCLS" --p-n-jobs "$THREADS" --p-reads-per-batch "$READS_PER_BATCH" --p-confidence "$TAXONOMY_CONFIDENCE" --p-read-orientation auto --o-classification "$UOUT"; outputs "ITS_UNITE_${UNITE_VERSION}" "$UOUT"
echo "Taxonomie terminee pour 16S et ITS."
