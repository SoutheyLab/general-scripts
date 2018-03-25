GATK_EXOME_INTERVALS="/scratch/vh83/sandbox/gene_amplicon_bed_metrics/Broad.human.exome.b37.interval_list"
BSTP_BED="/scratch/vh83/sandbox/gene_amplicon_bed_metrics/BRA-STRAP_621717_100.final.roverfile_g37.numsort.bed"
other_vep="/usr/local/vep/90/ensembl-vep/cache"
reference="/projects/vh83/reference/genomes/b37/bwa_0.7.12_index/human_g1k_v37_decoy.fasta"

BASENAME=${1}
METRICS_FOLDER=${2}

# Cat amplicon_metrics files
cat ${METRICS_FOLDER}*.amplicon-metrics.* > ${BASENAME}_amplicon_metrics_all

# Pass to python script to calculate averages (will append _averages to the output name)
python scripts/calc_amplicon_averages.py ${BASENAME}_amplicon_metrics_all \
    > ${BASENAME}_amplicon_metrics_all_averages

# Calculate pct-exonic and average depth for all amplicons using GATK intervals and python
grep -v "^@" ${GATK_EXOME_INTERVALS} | \
    bedtools coverage -a ${BSTP_BED} -b stdin | cut -f 9 | \
    paste <(python scripts/calc_amplicon_averages.py ${BASENAME}_amplicon_metrics_all) - > \
    ${BASENAME}_amplicon_metrics_pct_avg

# Remove header, cut columns and pass to VEP for annotation
cat ${BASENAME}_amplicon_metrics_pct_avg | cut -f 1-4 > ${BASENAME}_for_vep.tsv 
vep --cache --dir_cache ${other_vep} \
    --assembly GRCh37 --refseq --offline \
    --fasta ${reference} \
    --sift b --polyphen b --symbol --numbers --biotype \
    --total_length --hgvs \
    --force_overwrite --flag_pick --no_stats \
    --tab \
    -i ${BASENAME}_for_vep.tsv \
    -o ${BASENAME}_vep-anno.tsv

# Combine vep and other
python scripts/recombine_vep_avg.py \
    ${BASENAME}_vep-anno.tsv \
    ${BASENAME}_amplicon_metrics_pct_avg \
    ${BASENAME}_amplicon_info_complete.tsv
