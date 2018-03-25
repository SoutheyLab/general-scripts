import sys
import re


def consume_hd(fs):
    '''Consume lines starting with # given a file stream'''
    for line in fs:
        if line[:2] == '##':
            continue
        else:
            header = line[:-1]
            break
    return header


def store_vep(vepf):
    header = consume_hd(vepf)
    vep_dict = {'header': header}
    for line in vepf:
        splitl = line[:-1].split('\t')
        if splitl[17] != '1': continue  # Skipping non-vep picks
        d_key = re.search(r'.+:(\d+)-.+', splitl[1]).group(1)  # Keying by start
        if d_key in vep_dict:
            quit('Non-unique amplicon start key')
        else:
            vep_dict[d_key] = splitl
    return vep_dict


def combine_info(vep_dict, avgf, outputf):
    avg_hd = 'chrom\tstart\tend\tallele_for-vep\tforward\treverse\tavg_depth\tfraction_exonic'
    new_hd = avg_hd + '\t' + vep_dict['header'] + '\n' 
    outputf.write(new_hd)
    for line in avgf:
        splitl = line[:-1].split('\t')
        s_key = splitl[1]
        newl = splitl + vep_dict[s_key] 
        outputf.write('\t'.join(newl) + '\n')


def main():
    ipv = sys.argv[1]
    ipa = sys.argv[2]
    output_pth = sys.argv[3]

    with open(ipv, 'r') as vepf, open(ipa, 'r') as avgf, open(output_pth, 'w') as outputf:
        vep_dict = store_vep(vepf)
        combine_info(vep_dict, avgf, outputf)


if __name__ == '__main__':
    main()

