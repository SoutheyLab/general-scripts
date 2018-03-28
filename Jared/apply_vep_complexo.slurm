#!/bin/bash
#SBATCH --job-name=annotation
#SBATCH --account=vh83
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --partition=com
#SBATCH --mem=50000
#SBATCH --cpus-per-task=12

# Loading VEP
module load vep/90
# module load htslib/1.7

# Variables
# From the config file
reference="/projects/vh83/reference/genomes/b37/bwa_0.7.12_index/human_g1k_v37_decoy.fasta"
other_vep="/usr/local/vep/90/ensembl-vep/cache"
brcaexpath="/projects/vh83/reference/annotation_databases/BRCA-Exchange/BRCA-exchange_accessed-180118/BRCA-exchange_accessed-180118.sort.vcf.gz"
gnomadpath="/projects/vh83/reference/annotation_databases/gnomAD/gnomad.exomes.r2.0.2.sites/gnomad.exomes.r2.0.2.sites.vcf.gz"
revelpath="/projects/vh83/reference/annotation_databases/REVEL/REVEL-030616/revel_all_chromosomes.vcf.gz"
maxentscanpath="/projects/vh83/reference/annotation_databases/MaxEntScan/MaxEntScan_accessed-240118"
exacpath="/projects/vh83/reference/annotation_databases/ExAC/ExAC_nonTCGA.r0.3.1/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz"
dbnsfppath="/projects/vh83/reference/annotation_databases/dbNSFP/dbNSFPv2.9.3-VEP/dbNSFP-2.9.3.gz"
dbscsnvpath="/projects/vh83/reference/annotation_databases/dbscSNV/dbscSNV1.0-VEP/dbscSNV.txt.gz"
caddpath="/projects/vh83/reference/annotation_databases/CADD/CADD-v1.3/1000G_phase3.tsv.gz"

# Input
input_vcf="path/to/filtered_vcf.vcf"


# # Re-zipping and indexing with bgzip/tabix
# zcat ${input_vcf} | bgzip -c > ${input_vcf%.vcf.gz}.bgz.vcf.gz
# tabix -p vcf ${input_vcf%.vcf.gz}.bgz.vcf.gz

date
# Annotate with vep
vep --cache --dir_cache ${other_vep} \
    --assembly GRCh37 --refseq --offline \
    --fasta ${reference} \
    --sift b --polyphen b --symbol --numbers --biotype --total_length --hgvs --format vcf \
    --vcf --force_overwrite --flag_pick --no_stats \
    --custom ${brcaexpath},brcaex,vcf,exact,0,Clinical_significance_ENIGMA,Comment_on_clinical_significance_ENIGMA,Date_last_evaluated_ENIGMA,Pathogenicity_expert,HGVS_cDNA,HGVS_Protein,BIC_Nomenclature \
    --custom ${gnomadpath},gnomAD,vcf,exact,0,AF_NFE,AN_NFE \
    --custom ${revelpath},RVL,vcf,exact,0,REVEL_SCORE \
    --plugin MaxEntScan,${maxentscanpath} \
    --plugin ExAC,${exacpath},AC,AN \
    --plugin dbNSFP,${dbnsfppath},REVEL_score,REVEL_rankscore \
    --plugin dbscSNV,${dbscsnvpath} \
    --plugin CADD,${caddpath} \
    --fork 12 \
    -i ${input_vcf%.vcf.gz}.bgz.vcf.gz \
    -o ${input_vcf%.vcf.gz}.bgz.vep-anno.vcf
date
