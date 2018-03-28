#!/bin/bash
#SBATCH --job-name=annotation
#SBATCH --account=vh83
#SBATCH --time=03:00:00
#SBATCH --ntasks=1
#SBATCH --partition=com
#SBATCH --mem=16000
#SBATCH --cpus-per-task=1

# Sourcing environment and loading modules
source /projects/vh83/reference/GEMINI/gemini_env/bin/activate  # Environment for modified GEMINI
module load vep/90
module load vt/0.57
module load gatk/3.7
module load htslib/1.7

# Variables
CORES=6
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
CURR=${input_vcf}


# 1. Re-zipping and indexing with bgzip/tabix
zcat ${CURR} | bgzip -c > ${CURR%.vcf.gz}.bgz.vcf.gz
tabix -p vcf ${CURR%.vcf.gz}.bgz.vcf.gz
CURR=${CURR%.vcf.gz}.bgz.vcf.gz


# 2. Genotype filter with gatk
java -jar -Xmx8g /usr/local/gatk/3.7/executable/GenomeAnalysisTK.jar \
    -T VariantFiltration -R ${reference} \
    -V ${CURR} \
    -G_filter "g.isHet() == 1 && 1.0 * AB < 20.0" \
    -G_filterName "AltFreqLow" \
    -G_filter "DP < 10.0" \
    -G_filterName "LowDP" | \
    bgzip -c > ${CURR%.bgz.vcf.gz}.gt-filter.vcf.gz
CURR=${CURR%.bgz.vcf.gz}.gt-filter.vcf.gz


# 3. Decompose multiallelic sites with vt
vt decompose -s ${CURR} | vt normalize -r ${reference} - | bgzip -c > ${CURR%.vcf.gz}.decomp.norm.vcf.gz
CURR=${CURR%.vcf.gz}.decomp.norm.vcf.gz


# # # # The following steps require multiple cores
# # 4. Annotate with vep
# vep --cache --dir_cache ${other_vep} \
#     --assembly GRCh37 --refseq --offline \
#     --fasta ${reference} \
#     --sift b --polyphen b --symbol --numbers --biotype --total_length --hgvs --format vcf \
#     --vcf --force_overwrite --flag_pick --no_stats \
#     --custom ${brcaexpath},brcaex,vcf,exact,0,Clinical_significance_ENIGMA,Comment_on_clinical_significance_ENIGMA,Date_last_evaluated_ENIGMA,Pathogenicity_expert,HGVS_cDNA,HGVS_Protein,BIC_Nomenclature \
#     --custom ${gnomadpath},gnomAD,vcf,exact,0,AF_NFE,AN_NFE \
#     --custom ${revelpath},RVL,vcf,exact,0,REVEL_SCORE \
#     --plugin MaxEntScan,${maxentscanpath} \
#     --plugin ExAC,${exacpath},AC,AN \
#     --plugin dbNSFP,${dbnsfppath},REVEL_score,REVEL_rankscore \
#     --plugin dbscSNV,${dbscsnvpath} \
#     --plugin CADD,${caddpath} \
#     --fork ${CORES} \
#     -i ${CURR} \
#     -o ${CURR%.vcf.gz}.vep-anno.vcf


# # 5. Find and replace decomposed sites
# sed -E -e 's|([0-9])/\.|\1/0|g' -e 's|\./([0-9])|0/\1|g' ${CURR} | bgzip -c > ${CURR%.vcf}.f-and-r.vcf.gz
# tabix -p vcf ${CURR%.vcf}.f-and-r.vcf.gz
# CURR=${CURR%.vcf}.f-and-r.vcf.gz
# echo "Sed done."


# 6. Load the vcf with GEMINI (note: annotation blocks not decomposed, think it's a bad idea in this context)
gemini load --skip-gerp-bp --skip-cadd -t VEP --cores ${CORES} -v ${CURR} ${CURR%.vcf.gz}.db
echo "GEMINI done."
