#!/usr/bin/env bash
# =============================================================================
# 001_pipeline_QIIME2_PE_Pycnandra_V2.sh
# Projet : Pycnandra_V2
# Marqueurs : 16S et ITS, paired-end
#
# Selection explicite : seuls les echantillons biologiques *-5-* sont analyses.
# Les echantillons *-10-* sont exclus des le manifest et ne sont jamais importes.
# PYC-NEG est conserve jusqu'a DADA2 puis utilise dans 004 pour le retrait des
# ASV detectes dans les controles negatifs. Il n'entre pas dans les analyses
# de diversite une fois la table decontaminee produite.
# =============================================================================
set -Eeuo pipefail

export JAVA_HOME="${JAVA_HOME:-}"
export JAVA_LD_LIBRARY_PATH="${JAVA_LD_LIBRARY_PATH:-}"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

shopt -s nullglob
IFS=$'\n\t'
trap 'rc=$?; echo "[ERREUR] Code ${rc}, ligne ${LINENO}: ${BASH_COMMAND}" >&2; exit "${rc}"' ERR

PROJECT_NAME="Pycnandra_V2"
PROJECT_DIR="/nvme/bio/data_fungi/${PROJECT_NAME}"
RAW_ROOT_DIR="${PROJECT_DIR}/01_raw_data"
RESULTS_DIR="${PROJECT_DIR}/02_amplicon_pipeline"
TMPDIR_BASE="${PROJECT_DIR}/tmp"
LOG_DIR="${RESULTS_DIR}/logs"

THREADS=1
QIIME_THREADS=8
TRIMMOMATIC_HEAP="4G"
FASTQC_ENV="fastqc"
MULTIQC_ENV="multiqc"
TRIMMOMATIC_ENV="trimmomatic"
QIIME2_ENV="qiime2-amplicon-2025.7"
PYTHON_ENV="excel_tools"

TRIMMOMATIC_ADAPTERS=true
ADAPTER_FILE="${PROJECT_DIR}/99_softwares/adapters_sequences.fasta"
LEADING=30
TRAILING=30
SLIDINGWINDOW="26:30"
MINLEN=150

RUN_FASTQC_RAW=true
RUN_TRIMMOMATIC=true
RUN_FASTQC_CLEAN=true
RUN_QIIME_IMPORT=true
RUN_DADA2=true
RUN_TREE=true
RUN_EXPORT=true
MARKERS=("16S" "ITS")

# A ajuster apres inspection de demux.qzv. 0 signifie : aucune troncature QIIME2.
DADA2_TRIM_LEFT_F_16S=0
DADA2_TRIM_LEFT_R_16S=0
DADA2_TRUNC_LEN_F_16S=0
DADA2_TRUNC_LEN_R_16S=0
DADA2_MAX_EE_F_16S=2
DADA2_MAX_EE_R_16S=2
DADA2_CHIM_METHOD_16S="consensus"
DADA2_TRIM_LEFT_F_ITS=0
DADA2_TRIM_LEFT_R_ITS=0
DADA2_TRUNC_LEN_F_ITS=0
DADA2_TRUNC_LEN_R_ITS=0
DADA2_MAX_EE_F_ITS=2
DADA2_MAX_EE_R_ITS=2
DADA2_CHIM_METHOD_ITS="consensus"

log(){ printf '[%(%F %T)T] %s\n' -1 "$*" | tee -a "${LOG_DIR}/pipeline.log"; }
die(){ log "ERREUR : $*"; exit 1; }
activate_env(){
  export JAVA_HOME="${JAVA_HOME:-}"
  export JAVA_LD_LIBRARY_PATH="${JAVA_LD_LIBRARY_PATH:-}"
  export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
  conda deactivate >/dev/null 2>&1 || true
  conda activate "$1"
}
check_command(){ command -v "$1" >/dev/null 2>&1 || die "Commande introuvable : $1"; }
marker_to_lower(){ tr '[:upper:]' '[:lower:]' <<< "$1"; }

set_marker_variables(){
  local marker="$1" marker_lc
  marker_lc="$(marker_to_lower "$marker")"
  RAW_DIR="${RAW_ROOT_DIR}/${marker}"
  MARKER_DIR="${RESULTS_DIR}/${marker_lc}"
  QC_RAW_DIR="${MARKER_DIR}/01_qc_raw"
  CLEAN_DIR="${MARKER_DIR}/02_cleaned_data"
  QC_CLEAN_DIR="${MARKER_DIR}/03_qc_cleaned"
  DATABASE_DIR="${MARKER_DIR}/04_database_files"
  QIIME_DIR="${MARKER_DIR}/05_qiime2"
  QIIME_CORE="${QIIME_DIR}/core"
  QIIME_VISUAL="${QIIME_DIR}/visual"
  QIIME_TREE="${QIIME_DIR}/tree"
  QIIME_EXPORT="${QIIME_DIR}/export"
  SAMPLE_SHEET="${DATABASE_DIR}/samples_${marker_lc}.tsv"
  RAW_FASTQ_LIST="${DATABASE_DIR}/raw_fastq_list_${marker_lc}.txt"
  RAW_PAIRS_TSV="${DATABASE_DIR}/raw_pairs_${marker_lc}.tsv"
  MANIFEST="${DATABASE_DIR}/manifest_pe_${marker_lc}.tsv"
  METADATA="${DATABASE_DIR}/sample-metadata_${marker_lc}.tsv"
  case "$marker" in
    16S) DADA2_TRIM_LEFT_F="$DADA2_TRIM_LEFT_F_16S"; DADA2_TRIM_LEFT_R="$DADA2_TRIM_LEFT_R_16S"; DADA2_TRUNC_LEN_F="$DADA2_TRUNC_LEN_F_16S"; DADA2_TRUNC_LEN_R="$DADA2_TRUNC_LEN_R_16S"; DADA2_MAX_EE_F="$DADA2_MAX_EE_F_16S"; DADA2_MAX_EE_R="$DADA2_MAX_EE_R_16S"; DADA2_CHIM_METHOD="$DADA2_CHIM_METHOD_16S";;
    ITS) DADA2_TRIM_LEFT_F="$DADA2_TRIM_LEFT_F_ITS"; DADA2_TRIM_LEFT_R="$DADA2_TRIM_LEFT_R_ITS"; DADA2_TRUNC_LEN_F="$DADA2_TRUNC_LEN_F_ITS"; DADA2_TRUNC_LEN_R="$DADA2_TRUNC_LEN_R_ITS"; DADA2_MAX_EE_F="$DADA2_MAX_EE_F_ITS"; DADA2_MAX_EE_R="$DADA2_MAX_EE_R_ITS"; DADA2_CHIM_METHOD="$DADA2_CHIM_METHOD_ITS";;
  esac
}

run_fastqc_multiqc_from_list(){
  local list="$1" outdir="$2" label="$3"; local files=()
  mapfile -t files < "$list"; ((${#files[@]})) || die "Liste FastQC vide : $list"
  mkdir -p "$outdir/fastqc" "$outdir/multiqc"
  activate_env "$FASTQC_ENV"; check_command fastqc
  fastqc --threads "$THREADS" --outdir "$outdir/fastqc" "${files[@]}"
  activate_env "$MULTIQC_ENV"; check_command multiqc
  multiqc --force --outdir "$outdir/multiqc" "$outdir/fastqc"
  log "FastQC/MultiQC termine : $label"
}
run_fastqc_multiqc_from_directory(){
  local indir="$1" outdir="$2" label="$3"; local files=("$indir"/*_paired.fastq.gz)
  ((${#files[@]})) || die "Aucun FASTQ paired dans $indir"
  mkdir -p "$outdir/fastqc" "$outdir/multiqc"
  activate_env "$FASTQC_ENV"; check_command fastqc
  fastqc --threads "$THREADS" --outdir "$outdir/fastqc" "${files[@]}"
  activate_env "$MULTIQC_ENV"; check_command multiqc
  multiqc --force --outdir "$outdir/multiqc" "$outdir/fastqc"
  log "FastQC/MultiQC termine : $label"
}

CONDA_BASE="$(conda info --base 2>/dev/null || true)"
[[ -n "$CONDA_BASE" && -f "$CONDA_BASE/etc/profile.d/conda.sh" ]] || { echo "Conda introuvable" >&2; exit 1; }
source "$CONDA_BASE/etc/profile.d/conda.sh"
mkdir -p "$RESULTS_DIR" "$LOG_DIR" "$TMPDIR_BASE"; export TMPDIR="$TMPDIR_BASE"
[[ -d "$RAW_ROOT_DIR" ]] || die "Raw root absent : $RAW_ROOT_DIR"
log "Demarrage ${PROJECT_NAME}; seuls les echantillons -5- et PYC-NEG sont retenus."

# 1. Inventaire, selection stricte et creation manifest/metadata
activate_env "$PYTHON_ENV"; check_command python
for marker in "${MARKERS[@]}"; do
  set_marker_variables "$marker"
  mkdir -p "$QC_RAW_DIR" "$CLEAN_DIR" "$QC_CLEAN_DIR" "$DATABASE_DIR" "$QIIME_CORE" "$QIIME_VISUAL" "$QIIME_TREE" "$QIIME_EXPORT"
  python - "$RAW_DIR" "$marker" "$SAMPLE_SHEET" "$RAW_FASTQ_LIST" "$RAW_PAIRS_TSV" "$MANIFEST" "$METADATA" "$CLEAN_DIR" <<'PY'
import re, sys
from pathlib import Path
import pandas as pd
raw_dir, marker, sample_sheet, raw_list, raw_pairs_file, manifest, metadata, clean_dir = map(Path, sys.argv[1:])
pat=re.compile(r'^(?P<sample>.+)_S\d+_L\d{3}_R(?P<read>[12])_\d{3}\.fastq\.gz$')
records={}
for f in sorted(raw_dir.glob('*.fastq.gz')):
    m=pat.match(f.name)
    if not m: continue
    sample=m.group('sample')
    # Exclusion stricte de tous les -10- ; retention des biologiques -5- et PYC-NEG.
    if '-10-' in sample: continue
    if not ('-5-' in sample or sample == 'PYC-NEG'):
        continue
    if sample == 'PYC-NEG':
        kind='negative_control'; host='NA'; replicate='NA'; dilution='NA'
    else:
        x=re.fullmatch(r'(?P<host>Pa|Pf)-(?P<rep>\d+)-(?P<dilution>5)-(?P<library>\d+)', sample)
        if not x: raise SystemExit(f'Nom biologique inattendu: {sample}')
        kind='biological'; host=x['host']; replicate=x['rep']; dilution=x['dilution']
    records.setdefault(sample, {'sample-id':sample,'marker':marker,'sample_type':kind,'host_code':host,'replicate':replicate,'dilution':dilution})[f'R{m["read"]}']=str(f.resolve())
rows=[]
for sample,d in sorted(records.items()):
    if 'R1' not in d or 'R2' not in d: raise SystemExit(f'Paire incomplete: {sample}')
    rows.append(d)
df=pd.DataFrame(rows)
if df.empty: raise SystemExit(f'Aucun echantillon retenu dans {raw_dir}')
if set(df['sample-id']) != {'PYC-NEG','Pa-1-5-1','Pa-2-5-1','Pa-3-5-1','Pa-4-5-1','Pa-5-5-1','Pa-6-5-1','Pf-1-5-1','Pf-2-5-1','Pf-3-5-1','Pf-4-5-1','Pf-5-5-1','Pf-6-5-1'}:
    raise SystemExit('Jeu retenu inattendu: verifier les noms ou la selection.')
df.to_csv(sample_sheet,sep='\t',index=False)
pairs=df[['sample-id','R1','R2']].rename(columns={'R1':'raw-forward-absolute-filepath','R2':'raw-reverse-absolute-filepath'})
pairs.to_csv(raw_pairs_file,sep='\t',index=False)
with open(raw_list,'w') as h:
    for p in list(pairs.iloc[:,1])+list(pairs.iloc[:,2]): h.write(p+'\n')
pd.DataFrame({'sample-id':df['sample-id'],'forward-absolute-filepath':[str(clean_dir/f'{s}_R1_paired.fastq.gz') for s in df['sample-id']],'reverse-absolute-filepath':[str(clean_dir/f'{s}_R2_paired.fastq.gz') for s in df['sample-id']]}).to_csv(manifest,sep='\t',index=False)
df.drop(columns=['R1','R2']).rename(columns={'sample-id':'#SampleID'}).to_csv(metadata,sep='\t',index=False)
print(f'{marker}: {len(df)} echantillons retenus (12 biologiques -5- + PYC-NEG)')
PY
  log "$marker: manifest et metadata generes"
done

# 2. QC raw
if [[ "$RUN_FASTQC_RAW" == true ]]; then for marker in "${MARKERS[@]}"; do set_marker_variables "$marker"; run_fastqc_multiqc_from_list "$RAW_FASTQ_LIST" "$QC_RAW_DIR" "raw $marker"; done; fi

# 3. Trimmomatic
if [[ "$RUN_TRIMMOMATIC" == true ]]; then
  [[ "$TRIMMOMATIC_ADAPTERS" == false || -f "$ADAPTER_FILE" ]] || die "Adaptateurs absents: $ADAPTER_FILE"
  activate_env "$TRIMMOMATIC_ENV"; check_command trimmomatic
  for marker in "${MARKERS[@]}"; do
    set_marker_variables "$marker"
    while IFS=$'\t' read -r id r1 r2; do
      [[ "$id" == sample-id ]] && continue
      o1="$CLEAN_DIR/${id}_R1_paired.fastq.gz"; u1="$CLEAN_DIR/${id}_R1_unpaired.fastq.gz"; o2="$CLEAN_DIR/${id}_R2_paired.fastq.gz"; u2="$CLEAN_DIR/${id}_R2_unpaired.fastq.gz"
      if [[ -s "$o1" && -s "$o2" ]]; then continue; fi
      args=(); [[ "$TRIMMOMATIC_ADAPTERS" == true ]] && args=("ILLUMINACLIP:${ADAPTER_FILE}:2:30:10")
      trimmomatic PE -Xmx"$TRIMMOMATIC_HEAP" -threads "$THREADS" -phred33 "$r1" "$r2" "$o1" "$u1" "$o2" "$u2" "${args[@]}" "LEADING:${LEADING}" "TRAILING:${TRAILING}" "SLIDINGWINDOW:${SLIDINGWINDOW}" "MINLEN:${MINLEN}"
    done < "$RAW_PAIRS_TSV"
  done
fi

# 4. QC cleaned
if [[ "$RUN_FASTQC_CLEAN" == true ]]; then for marker in "${MARKERS[@]}"; do set_marker_variables "$marker"; run_fastqc_multiqc_from_directory "$CLEAN_DIR" "$QC_CLEAN_DIR" "cleaned $marker"; done; fi

# 5. QIIME2 : import, DADA2, arbre, exports
activate_env "$QIIME2_ENV"; check_command qiime
for marker in "${MARKERS[@]}"; do
  set_marker_variables "$marker"
  [[ -s "$MANIFEST" && -s "$METADATA" ]] || die "Manifest/metadata absents pour $marker"
  if [[ "$RUN_QIIME_IMPORT" == true ]]; then
    rm -f "$QIIME_CORE/demux.qza" "$QIIME_VISUAL/demux.qzv"
    qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path "$MANIFEST" --input-format PairedEndFastqManifestPhred33V2 --output-path "$QIIME_CORE/demux.qza"
    qiime demux summarize --i-data "$QIIME_CORE/demux.qza" --o-visualization "$QIIME_VISUAL/demux.qzv"
  fi
  if [[ "$RUN_DADA2" == true ]]; then
    rm -f "$QIIME_CORE/table.qza" "$QIIME_CORE/rep-seqs.qza" "$QIIME_CORE/denoising-stats.qza"
    qiime dada2 denoise-paired --i-demultiplexed-seqs "$QIIME_CORE/demux.qza" --p-trim-left-f "$DADA2_TRIM_LEFT_F" --p-trim-left-r "$DADA2_TRIM_LEFT_R" --p-trunc-len-f "$DADA2_TRUNC_LEN_F" --p-trunc-len-r "$DADA2_TRUNC_LEN_R" --p-max-ee-f "$DADA2_MAX_EE_F" --p-max-ee-r "$DADA2_MAX_EE_R" --p-chimera-method "$DADA2_CHIM_METHOD" --p-n-threads "$QIIME_THREADS" --o-table "$QIIME_CORE/table.qza" --o-representative-sequences "$QIIME_CORE/rep-seqs.qza" --o-denoising-stats "$QIIME_CORE/denoising-stats.qza"
    qiime metadata tabulate --m-input-file "$QIIME_CORE/denoising-stats.qza" --o-visualization "$QIIME_VISUAL/denoising-stats.qzv"
    qiime feature-table summarize --i-table "$QIIME_CORE/table.qza" --m-sample-metadata-file "$METADATA" --o-visualization "$QIIME_VISUAL/table.qzv"
    qiime feature-table tabulate-seqs --i-data "$QIIME_CORE/rep-seqs.qza" --o-visualization "$QIIME_VISUAL/rep-seqs.qzv"
  fi
  if [[ "$RUN_TREE" == true ]]; then
    qiime alignment mafft --i-sequences "$QIIME_CORE/rep-seqs.qza" --p-n-threads "$QIIME_THREADS" --o-alignment "$QIIME_TREE/aligned-rep-seqs.qza"
    qiime alignment mask --i-alignment "$QIIME_TREE/aligned-rep-seqs.qza" --o-masked-alignment "$QIIME_TREE/masked-aligned-rep-seqs.qza"
    qiime phylogeny fasttree --i-alignment "$QIIME_TREE/masked-aligned-rep-seqs.qza" --o-tree "$QIIME_TREE/unrooted-tree.qza"
    qiime phylogeny midpoint-root --i-tree "$QIIME_TREE/unrooted-tree.qza" --o-rooted-tree "$QIIME_TREE/rooted-tree.qza"
  fi
  if [[ "$RUN_EXPORT" == true ]]; then
    rm -rf "$QIIME_EXPORT/core" "$QIIME_EXPORT/tree" "$QIIME_EXPORT/visual"; mkdir -p "$QIIME_EXPORT/core" "$QIIME_EXPORT/tree" "$QIIME_EXPORT/visual"
    qiime tools export --input-path "$QIIME_CORE/table.qza" --output-path "$QIIME_EXPORT/core/table"
    qiime tools export --input-path "$QIIME_CORE/rep-seqs.qza" --output-path "$QIIME_EXPORT/core/rep-seqs"
    qiime tools export --input-path "$QIIME_CORE/denoising-stats.qza" --output-path "$QIIME_EXPORT/core/denoising-stats"
    qiime tools export --input-path "$QIIME_TREE/rooted-tree.qza" --output-path "$QIIME_EXPORT/tree/rooted-tree"
  fi
  log "$marker termine"
done
log "Pipeline 001 termine. Lancer ensuite 003 puis 002, puis 004."
