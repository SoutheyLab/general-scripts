# Splits variants into multiple lines (one for each VEP annotation block, if there are multiple)
import sys
import re

input_vcf_path = sys.argv[1]
output_decomposed_vcf_path = input_vcf_path[:-4] + ".anno-decomp.vcf"

with open(input_vcf_path, 'r') as input_vcf, open(output_decomposed_vcf_path, 'w') as output_vcf:
    for line in input_vcf:
        if line.startswith('#'):
            output_vcf.write(line)
            continue
        splitline = line.rstrip('\n').split('\t')
        infofield = splitline[7]
        if "CSQ=" in infofield:
            vep_block = re.search("CSQ=(.+)(?:$|;)", infofield).group(1)
        else:
            output_vcf.write(line)
            print "No VEP info block for {chrom}:{pos}".format(chrom=splitline[0], pos=splitline[1])
            continue
        if ',' in vep_block:
            info_strip_vep = [field for field in infofield.split(';') if not field.startswith("CSQ=")]
            for annoblock in vep_block.split(','):
                new_info_field = info_strip_vep + ["CSQ=" + annoblock]
                finished_line = splitline[:7] + [';'.join(new_info_field)] + splitline[8:]
                output_vcf.write('\t'.join(finished_line) + '\n')
        else:
            output_vcf.write(line)
