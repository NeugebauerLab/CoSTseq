import pysam
import numpy as np
from src.bam_utils import read_pair_generator, read_ends
import matplotlib.pyplot as plt
import argparse
from os.path import exists
plt.rcParams['pdf.fonttype'] = 42

def insert_length(files):
    
    for n, filename in enumerate(files):

        assert exists(filename + '.bai'), "Index file not in directory."

        bamfile = pysam.AlignmentFile(filename, 'rb')
        insert_lengths = []
        
        try:

            for read1, read2 in read_pair_generator(bamfile):

                _, fivep_r1, _ = read_ends(read1)
                _, fivep_r2, _ = read_ends(read2)

                insert_lengths.append(abs(fivep_r1-fivep_r2))

            fig = plt.figure()
            plt.hist(insert_lengths, bins=np.arange(0, 1200, 10))
            plt.title(f"mean: {np.mean(insert_lengths):.2f}, median: {np.median(insert_lengths):.2f}, n: {len(insert_lengths)}")
            plt.xlabel('insert size')
            plt.ylabel('number of reads')
            plt.ioff()
            fig.savefig(filename[:-3] + 'pdf')

            print(f"Done with file {n+1}.")

        finally:
            bamfile.close()



if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('bam', nargs='+', help='The paired-end BAM file that needs to be analyzed.')
    args = parser.parse_args()

    insert_length(args.bam)

