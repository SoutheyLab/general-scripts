import pysam
import sys


def main():
    # bamsuffix = '.sorted.locatit.bam'
    bam_f = pysam.AlignmentFile(sys.argv[1], 'rb')
    output_pth = sys.argv[2]

    counts = []

    with open(output_pth, 'w') as output_f:
        for line in bam_f.fetch():
            counts.append(str(line.get_tag('XI')))
        output_f.write('\n'.join(counts))



if __name__ == '__main__':
    main()
