import pysam
import pandas as pd
import numpy as np
from src.bam_utils import read_pair_generator
import re
from tqdm import tqdm
import argparse
import os.path
import tempfile
from shutil import which
import subprocess

def deduplicate(file, UMI_prefix, UMI_len, oname, without_umi, ignore_seq):

    # make sure samtools exists
    assert which('samtools') is not None, "samtools not installed or not in PATH"

    fd, path = tempfile.mkstemp()

    # open bam file
    bam = pysam.AlignmentFile(file, 'rb')
    chrom_names = bam.header.references

    try: # decorators for writing to temp file
        with os.fdopen(fd, 'w') as out:
            
            # for statistics
            n_rem = 0
            n_tot = 0
            
            # iterate over chromosomes to save memory
            # for larger chromosomes or more reads this still takes a lot of memory, could consider splitting up the chromosomes
            for chrom in tqdm(chrom_names):

                UMI = []
                STT_R1 = []
                STT_R2 = []
                SEQ_R1 = []
                SEQ_R2 = []
                RID_R1 = []
                RID_R2 = []

                # append meaningful attributes of each read pair to a list
                for read1, read2 in read_pair_generator(bam, reference=chrom):
                    # append UMI if exists, otherwise this column will be non-informative
                    if without_umi:
                        UMI.append(0)
                    else:
                        UMI.append(re.findall(f'{UMI_prefix}[ACGTN]{{{UMI_len}}}', read1.query_name)[0].strip(UMI_prefix))
                    
                    STT_R1.append(read1.reference_start)
                    STT_R2.append(read2.reference_start)
                    SEQ_R1.append(read1.get_forward_sequence().upper())
                    SEQ_R2.append(read2.get_forward_sequence().upper())
                    RID_R1.append(read1.query_name)
                    RID_R2.append(read2.query_name)

                # detect duplicates based on relevant attributes
                comp_R1 = pd.DataFrame({'RID': RID_R1, 'UMI': UMI, 'STT': STT_R1, 'SEQ': SEQ_R1})
                comp_R2 = pd.DataFrame({'RID': RID_R2, 'UMI': UMI, 'STT': STT_R2, 'SEQ': SEQ_R2})

                if not ignore_seq:
                    dup_mask = comp_R1.duplicated(['UMI', 'STT', 'SEQ'], keep='first') & comp_R2.duplicated(['UMI', 'STT', 'SEQ'], keep='first')
                else:
                    dup_mask = comp_R1.duplicated(['UMI', 'STT'], keep='first') & comp_R2.duplicated(['UMI', 'STT'], keep='first')

                # logging
                n_rem += np.sum(dup_mask)
                n_tot += len(dup_mask)

                # write a list of read names to keep, append to it for each chromosome
                # would be better to write the reads that need to be removed (smaller file), but samtools doesnt't have
                # an option for removing reads based on a list.
                keep = list(comp_R1[~dup_mask]['RID']) + list(comp_R2[~dup_mask]['RID'])            
                for read in keep:
                    out.write(read + '\n')
                    
        bam.close()

        if n_rem != 0:

            print(f'Found {n_rem} duplicate reads ({100*n_rem/n_tot:.2f}%).')
            print('Filtering bam...')

            out = subprocess.run(['samtools', 'view', '-b',
                                  '-N', path,
                                  '-o', oname,
                                  file],
                                  capture_output=False, shell=False)

            assert out.returncode==0, "samtools view run failed."
            print('Done.')

        else:
            print('No duplicates found. Output file is identical to input.')

            out = subprocess.run(['cp', file, oname], capture_output=False, shell=False)


    finally: # close temp file
        os.remove(path)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('bam', help='The paired-end BAM file that needs to be de-duplicated.')
    parser.add_argument('-umi', help='UMI prefix that was appended to the read name (e.g. by fastp). Default is "UMI_"')
    parser.add_argument('-umi_len', help='Length of UMI. Default is 5.', type=int)
    parser.add_argument('-o', help='Output file name, optional.')
    parser.add_argument('--without_umi', help='Deduplicate only on sequence identity and alignment position.', action='store_true')
    parser.add_argument('--ignore_seq', help='Deduplicate only on UMI and alignment position.', action='store_true')
    args = parser.parse_args()

    # handle input args
    if not args.umi:
        UMI_prefix = 'UMI_'
    else:
        UMI_prefix = args.umi

    if not args.umi_len:
        UMI_len = 5
    else:
        UMI_len = args.umi_len

    # handle output file
    if not args.o:
        oname = file[:-4] + '_dedup.bam'
    else:
        oname = args.o

    deduplicate(args.bam, UMI_prefix, UMI_len, oname, args.without_umi, args.ignore_seq)

