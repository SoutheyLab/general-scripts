#!/bin/bash
#SBATCH --job-name=gemini
#SBATCH --account=vh83
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --partition=short
#SBATCH --mem=32000
#SBATCH --cpus-per-task=6

# Script to preprocess and load a VEP annotated vcf into a GEMINI database
# To run: sbatch load_mod_gemini.slurm --export=INVCF="/path/to/input/vcf.vcf",ALL

# Module loads
module load htslib/1.7

# Variables
DECOMP="/home/jburke/scripts/github/SoutheyLab/general-scripts/Jared/decompose_annotation_blocks.py"
CURRINPUT=${INVCF}

# Source the environment for the modified GEMINI (to read in FT information)
source /projects/vh83/reference/GEMINI/gemini_env/bin/activate


# Decompose annotation blocks (optional) (outputs with suffix .anno-decomp.vcf)
python ${DECOMP} ${CURRINPUT}
CURRINPUT=${CURRINPUT%.vcf}.anno-decomp.vcf


# Fix decomposed multiallelic sites
sed -E -e 's|([0-9])/\.|\1/0|g' -e 's|\./([0-9])|0/\1|g' ${CURRINPUT} > \
    ${CURRINPUT%.vcf}.f-and-r.vcf
CURRINPUT=${CURRINPUT%.vcf}.f-and-r.vcf


# Load the database
gemini load --skip-gerp-bp --skip-cadd -t VEP --cores 6 -v ${CURRINPUT} ${CURRINPUT%.vcf}.db
