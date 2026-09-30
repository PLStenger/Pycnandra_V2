#!/usr/bin/env bash
# 002_rarefaction_Pycnandra_V2.sh
# Soustraction des ASV de PYC-NEG, puis courbes de rarefaction
# pour les seuls échantillons biologiques *-5-*.

set -Eeuo pipefail

PROJECT_DIR="/nvme/bio/data_fungi/Pycnandra_V2"
RESULTS_DIR="${PROJECT_DIR}/02_amplicon_pipeline"
QIIME_ENV="qiime2-amplicon-2025.7"
BIOM_ENV="biom-format"

eval "$(conda shell.bash hook)"

qiime_run() {
    conda run -n "${QIIME_ENV}" qiime "$@"
}

biom_run() {
    conda run -n "${BIOM_ENV}" biom "$@"
}

for marker in 16S ITS; do
    lc="$(tr '[:upper:]' '[:lower:]' <<< "${marker}")"

    QIIME_DIR="${RESULTS_DIR}/${lc}/05_qiime2"
    DATABASE="${RESULTS_DIR}/${lc}/04_database_files"
    DECON="${QIIME_DIR}/decontam"
    VISUAL="${QIIME_DIR}/visual"

    INPUT="${QIIME_DIR}/core/table.qza"
    TREE="${QIIME_DIR}/tree/rooted-tree.qza"
    META="${DATABASE}/sample-metadata_${lc}.tsv"

    NEG_TABLE="${DECON}/table_PYC_NEG_only.qza"
    NEG_EXPORT="${DECON}/export_PYC_NEG"
    NEG_TSV="${DECON}/table_PYC_NEG.tsv"
    NEG_IDS="${DECON}/negative_asv_ids.tsv"

    NO_NEG_ASV="${DECON}/table_no_negative_asvs.qza"
    BIO="${DECON}/table_no_negative_asvs_biological_only.qza"

    BIO_EXPORT="${DECON}/export_biological"
    BIO_TSV="${DECON}/table_biological_decontaminated.tsv"
    OUT="${VISUAL}/alpha-rarefaction_${lc}_decontam_biological_only.qzv"

    mkdir -p "${DECON}" "${VISUAL}"

    for file in "${INPUT}" "${TREE}" "${META}"; do
        [[ -s "${file}" ]] || {
            echo "${marker} : fichier manquant ou vide : ${file}" >&2
            exit 1
        }
    done

    echo "=== ${marker} : isolement de PYC-NEG ==="

    rm -f "${NEG_TABLE}" "${NO_NEG_ASV}" "${BIO}" "${OUT}"
    rm -rf "${NEG_EXPORT}" "${BIO_EXPORT}"

    qiime_run feature-table filter-samples \
        --i-table "${INPUT}" \
        --m-metadata-file "${META}" \
        --p-where "[sample_type]='negative_control'" \
        --o-filtered-table "${NEG_TABLE}"

    qiime_run tools export \
        --input-path "${NEG_TABLE}" \
        --output-path "${NEG_EXPORT}"

    biom_run convert \
        -i "${NEG_EXPORT}/feature-table.biom" \
        -o "${NEG_TSV}" \
        --to-tsv

       # Fichier de métadonnées QIIME2 : en-tête reconnu + IDs des ASV.
    {
        printf 'feature-id\n'
        awk -F '\t' '
            $1 !~ /^#/ &&
            length($1) == 32 &&
            $1 ~ /^[[:xdigit:]]+$/ {
                print $1
            }
        ' "${NEG_TSV}"
    } > "${NEG_IDS}"
    

    n_asv="$(awk 'END {print NR - 1}' "${NEG_IDS}")"
    (( n_asv > 0 )) || {
        echo "${marker} : aucun ASV détecté dans PYC-NEG ; arrêt pour vérification." >&2
        exit 1
    }

    echo "${marker} : ${n_asv} ASV détectés dans PYC-NEG."

    echo "=== ${marker} : soustraction des ASV du contrôle ==="

    qiime_run feature-table filter-features \
        --i-table "${INPUT}" \
        --m-metadata-file "${NEG_IDS}" \
        --p-exclude-ids \
        --o-filtered-table "${NO_NEG_ASV}"

    echo "=== ${marker} : exclusion de PYC-NEG de la table biologique ==="

    qiime_run feature-table filter-samples \
        --i-table "${NO_NEG_ASV}" \
        --m-metadata-file "${META}" \
        --p-where "[sample_type]='biological'" \
        --o-filtered-table "${BIO}"

    # Export de la table biologique pour obtenir les profondeurs réelles
    # après décontamination, sans reprendre celles du projet Araucaria.
    qiime_run tools export \
        --input-path "${BIO}" \
        --output-path "${BIO_EXPORT}"

    biom_run convert \
        -i "${BIO_EXPORT}/feature-table.biom" \
        -o "${BIO_TSV}" \
        --to-tsv

    # Somme des lectures par colonne (échantillon), puis maximum observé.
       max_depth="$(
        awk -F '\t' '
            $1 ~ /^#OTU ID$/ || $1 ~ /^#OTU ID[[:space:]]*$/ {
                n = NF
                next
            }
            $1 ~ /^#/ || NF < 2 {
                next
            }
            {
                for (i = 2; i <= NF; i++) {
                    total[i] += $i
                }
            }
            END {
                max = 0
                for (i = 2; i <= n; i++) {
                    if (total[i] > max) max = total[i]
                }
                printf "%.0f\n", max
            }
        ' "${BIO_TSV}"
    )"

    (( max_depth > 0 )) || {
        echo "${marker} : profondeur maximale nulle après décontamination." >&2
        exit 1
    }

    echo "${marker} : courbe de rarefaction, profondeur maximale = ${max_depth}"

    qiime_run diversity alpha-rarefaction \
        --i-table "${BIO}" \
        --i-phylogeny "${TREE}" \
        --p-min-depth 1 \
        --p-max-depth "${max_depth}" \
        --m-metadata-file "${META}" \
        --o-visualization "${OUT}"

    echo "${marker} : table décontaminée : ${BIO}"
    echo "${marker} : courbe de rarefaction : ${OUT}"
done

echo "002 terminé pour 16S et ITS."
