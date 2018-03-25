#!/bin/bash
#SBATCH --job-name=gemini
#SBATCH --account=vh83
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --partition=short
#SBATCH --mem=64000
#SBATCH --cpus-per-task=12

# Script to preprocess and load a VEP annotated vcf into a GEMINI database
# To run modify INVCF to the desired input vcf

# Setting up environment
source /projects/vh83/reference/GEMINI/gemini_env/bin/activate  # Environment for modified GEMINI
module load htslib/1.7

# Variables
DECOMP="/home/jburke/scripts/github/SoutheyLab/general-scripts/Jared/decompose_annotation_blocks.py"
INVCF="/path/to/input/vcf.vcf"
CURRINPUT=${INVCF}


# Decompose annotation blocks (optional) (outputs with suffix .anno-decomp.vcf)
python ${DECOMP} ${CURRINPUT}
CURRINPUT=${CURRINPUT%.vcf}.anno-decomp.vcf
echo "Decomposition done."


# Fix decomposed multiallelic sites
sed -E -e 's|([0-9])/\.|\1/0|g' -e 's|\./([0-9])|0/\1|g' ${CURRINPUT} > \
    ${CURRINPUT%.vcf}.f-and-r.vcf
CURRINPUT=${CURRINPUT%.vcf}.f-and-r.vcf
echo "Sed done."


# Load the database
gemini load --skip-gerp-bp --skip-cadd -t VEP --cores 12 -v ${CURRINPUT} ${CURRINPUT%.vcf}.db
echo "GEMINI done."
