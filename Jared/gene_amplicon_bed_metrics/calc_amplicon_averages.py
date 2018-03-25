import sys


def calculate_averages(inputf, outputf):
    line_dict = {}
    output_lines = []
    for line in inputf:
        splitl = line[:-1].split('\t')
        d_key = splitl[3]  # Keying by forward primer
        reads = splitl[5]  # Number of reads at this location for this sample
        try:
            line_dict[d_key].append(float(reads))
        except KeyError:
            line_dict[d_key] = [float(reads)]
        bareline = '\t'.join(splitl[:-5])  # Saving the bed lines for later
        if bareline not in output_lines:  # Inefficient but can't be bothered
            output_lines.append(bareline)

    # outputf.write('chrom\tchromStart\tchromEnd\tfeature\tforward\treverse\taverage_depth\n')
    for line in output_lines:
        splitl = line.split('\t')
        d_key = splitl[3]
        avg = [str(float(sum(line_dict[d_key]))/float(len(line_dict[d_key])))]
        # Pretending they're all structural variants so VEP will annotate them, cheating 
        # but should work
        outline = splitl[:3] + ['DEL'] + splitl[3:] + avg
        outputf.write('\t'.join(outline) + '\n')


def main():
    input_pth = sys.argv[1]

    with open(input_pth, 'r') as inputf:
        calculate_averages(inputf, sys.stdout)


if __name__ == '__main__':
    main()
