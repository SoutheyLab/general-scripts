import gzip
import sys
import time
import numpy
import pstats, cProfile


def iterate_vcf(vcfin, vcfout):
    """Iterates through vcf."""
    var_dict = {}
    lc = 0
    smp2idx = {}
    curridx = 0
    hd_lines = {}
    # Consume and store the vcf header
    for line in vcfin:
        if line.startswith('#CHROM'):
            hd_lines['main'] = line
            break
        elif line.startswith('##FILTER'):
            if 'FILTER' in hd_lines:
                hd_lines['FILTER'].append(line)
            else:
                hd_lines['FILTER'] = [line]
        elif line.startswith('##INFO'):
            if 'INFO' in hd_lines:
                hd_lines['INFO'].append(line)
            else:
                hd_lines['INFO'] = [line]
        elif line.startswith('##contig'):
            if 'contig' in hd_lines:
                hd_lines['contig'].append(line)
            else:
                hd_lines['contig'] = [line]

        # These should be a single line each
        elif line.startswith('##reference'):
            hd_lines['reference'] = line
        elif line.startswith('##source'):
            hd_lines['source'] = line
        elif line.startswith('##fileformat'):
            hd_lines['fileformat'] = line
        elif line.startswith('##fileDate'):
            hd_lines['fileDate'] = line

        # Everything else should be commands or something
        else:
            if 'other' in hd_lines:
                hd_lines['other'].append(line)
            else:
                hd_lines['other'] = [line]
    # Fill dictionary
    tabjoin = '\t'.join
    for line in vcfin:
        # Splitting line
        splitl = line[:-1].split('\t')
        # if splitl[6][0] != 'P': continue  # Skip non-pass variants
        varkey = tabjoin(splitl[:5])
        inf = splitl[7].split(';')
        smp = inf[0][7:]
        # Referencing sample by index to save memory
        try:
            smpidx = smp2idx[smp]
        except KeyError:
            smp2idx[smp] = curridx
            smpidx = curridx
            curridx += 1
        try:
            var_dict[varkey].append(tabjoin([str(smpidx), inf[1][3:], inf[2][3:], inf[3][4:]]))
        except KeyError:
            var_dict[varkey] = [tabjoin([str(smpidx), inf[1][3:], inf[2][3:], inf[3][4:]])]
        # # Break conditions
        # lc += 1
        # if lc == 1000000: break
    # print(len(smp2idx.keys()))
    return var_dict, smp2idx, hd_lines


def write_header(vcfout, hd_lines):
    '''Writes everything in the header up to the final line (with CHROM and samples etc.)
    Careful as it expects most of this to be present.'''
    # Top lines
    vcfout.write(hd_lines['fileformat'])
    vcfout.write(hd_lines['fileDate'])
    vcfout.write(hd_lines['source'])
    # INFO, FILTER and FORMAT
    vcfout.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n')
    vcfout.write('##FORMAT=<ID=AB,Number=1,Type=Float,Description="Allele balance for each het genotype">\n') 
    vcfout.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">\n')
    vcfout.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">\n')
    # vcfout.write('##FORMAT=<ID=FT,Number=.,Type=String,Description="Genotype-level filter">\n')
    vcfout.write('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">\n')
    vcfout.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    vcfout.write('##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">\n')
    # Contigs, reference and other
    comp_loop = [vcfout.write(line) for line in hd_lines['contig']]
    vcfout.write(hd_lines['reference'])
    comp_loop = [vcfout.write(line) for line in hd_lines['other']]


def write_vcf(vcfout, var_dict, hd_lines, smp2idx):
    # Creating a vcf-style header and inverting the dictionary (will need to update the 
    # header writing above depending on what's done here
    idx2smp = {smp2idx[smp]: smp for smp in smp2idx}
    vcfout.write("{par}\tFORMAT\t{smp}\n".format(par=hd_lines['main'].rstrip('\n'),
                                                 smp='\t'.join([idx2smp[n] for n in range(0, len(idx2smp.keys()))])))

    # Base FORMAT field from GATK vcfs GT:AB:AD:DP:FT:GQ:PL (example
    # gt_field_novar = "./.:.:0,0:0:PASS:.:0,0,0"
    gt_field_novar = "./.:.:0,0:0:.:0,0,0"  # No FT
    numsmp = len(smp2idx.keys())
    for var in var_dict:
        # Writing the first bit of the line
        # varline = "{varkey}\t.\t.\tDP=50\tGT:AB:AD:DP:FT:GQ:PL\t".format(varkey=var)
        varline = "{varkey}\t.\t.\tDP=50\tGT:AB:AD:DP:GQ:PL\t".format(varkey=var)  # No FT

        # Creating and writing the genotype info
        gt_arr = numpy.empty(numsmp, dtype=object)
        gt_arr.fill(gt_field_novar)
        for smp in var_dict[var]:
            split_smp = smp.split('\t')
            idx = int(split_smp[0])
            nv = split_smp[1]
            np = split_smp[2]
            nr = str(int(np) - int(nv))
            pct = split_smp[3]
            # gt_field = "0/1:{ab}:{ref},{alt}:{dp}:PASS:99:0,0,0".format(ab=pct,
            #                                                             ref=nr,
            #                                                             alt=nv,
            #                                                             dp=np)
            gt_field = "0/1:{ab}:{ref},{alt}:{dp}:99:0,0,0".format(ab=pct,  # No FT
                                                                   ref=nr,
                                                                   alt=nv,
                                                                   dp=np)
            gt_arr[idx] = gt_field
        vcfout.write(''.join([varline, '\t'.join(gt_arr), '\n']))


def main():
    output_pth = sys.argv[1]
    with open(output_pth, 'wt') as outputf:
        var_dict, smp2idx, hd_lines = iterate_vcf(sys.stdin, outputf)
        write_header(outputf, hd_lines)
        write_vcf(outputf, var_dict, hd_lines, smp2idx)

    # # Profiling
    # cProfile.runctx("iterate_vcf(sys.stdin, None)", globals(), locals(), "Profile.prof")
    # s = pstats.Stats("Profile.prof")
    # s.strip_dirs().sort_stats("time").print_stats()

if __name__ == '__main__':
    main()

