#!/bin/bash
#SBATCH --job-name=annotation
#SBATCH --account=vh83
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --partition=short
#SBATCH --mem=110000
#SBATCH --cpus-per-task=24

# Loading VEP
module load vep/90
module load gatk/3.7
module load bcftools/1.7

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

# UNDR-ROVER / Input vcf
undr_rover_vcf="combined_undr_rover.vcf.gz"

# Run python script to convert the UNDR ROVER vcf to a standard format (all genotypes set as het)
# zcat ${undr_rover_vcf} | python collapse_reformat_undr_rover.py ${undr_rover_vcf%.vcf.gz}.pass.reformat.vcf

# bcftools sort -o ${undr_rover_vcf%.vcf.gz}.pass.reformat.sorted.vcf \
#     ${undr_rover_vcf%.vcf.gz}.pass.reformat.vcf

# java -jar -Xmx8g /usr/local/gatk/3.7/executable/GenomeAnalysisTK.jar \
#     -T VariantFiltration -R ${reference} \
#     -V ${undr_rover_vcf%.vcf.gz}.pass.reformat.sorted.vcf \
#     -o ${undr_rover_vcf%.vcf.gz}.pass.reformat.sorted.gt-filter.vcf \
#     -G_filter "g.isHet() == 1 && 1.0 * AB < 25.0" \
#     -G_filterName "AltFreqLow" \
#     -G_filter "DP < 25.0" \
#     -G_filterName "LowDP"


vep --cache --dir_cache ${other_vep} \
    --assembly GRCh37 --refseq --offline \
    --fasta ${reference} \
    --sift b --polyphen b --symbol --numbers --biotype --total_length --hgvs --format vcf \
    --vcf --force_overwrite --flag_pick --no_stats \
    --custom ${undr_rover_vcf},undrrover,vcf,exact,0,Sample,PCT,NV,NP \
    --custom ${brcaexpath},brcaex,vcf,exact,0,Clinical_significance_ENIGMA,Comment_on_clinical_significance_ENIGMA,Date_last_evaluated_ENIGMA,Pathogenicity_expert,HGVS_cDNA,HGVS_Protein,BIC_Nomenclature \
    --custom ${gnomadpath},gnomAD,vcf,exact,0,AF_NFE,AN_NFE \
    --custom ${revelpath},RVL,vcf,exact,0,REVEL_SCORE \
    --plugin MaxEntScan,${maxentscanpath} \
    --plugin ExAC,${exacpath},AC,AN \
    --plugin dbNSFP,${dbnsfppath},REVEL_score,REVEL_rankscore \
    --plugin dbscSNV,${dbscsnvpath} \
    --plugin CADD,${caddpath} \
    --fork 12 \
    -i ${undr_rover_vcf%.vcf.gz}.pass.reformat.sorted.gt-filter.vcf \
    -o ${undr_rover_vcf%.vcf.gz}.pass.reformat.sorted.gt-filter.vep-anno.vcf
