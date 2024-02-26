import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42
import pickle
from src.pro_utils import three_base_identity
from src.bam_utils import get_read_number
from src.DMS_Profile import DMS_Profile, PRO_Profile, Genome, Targeted_DMS_Profile
import json
from matplotlib.gridspec import GridSpec
import argparse

def QC(sample_name, output_path, fastp_json, mRNA_bam, mRNA_bam_ded, rRNA_bam, rRNA_bam_ded, dms_mRNA_pkl, dms_rRNA_pkl, genome_fasta, pro_mRNA_pkl):
    
    ## Initialize figure

    fig = plt.figure(layout="constrained", figsize=(10,6))

    gs = GridSpec(3, 3, figure=fig, height_ratios=(0.4, 0.3, 0.3))
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[0, 2], sharey=ax2)
    ax4 = fig.add_subplot(gs[1, 0])
    ax5 = fig.add_subplot(gs[2, 0], sharex=ax4, sharey=ax4)
    ax6 = fig.add_subplot(gs[1,1:])
    ax7 = fig.add_subplot(gs[2,1:], sharey=ax6)

    fig.suptitle(sample_name)

    ## read numbers-----------------------------------------------------------------------------------

    # grab total number of reads from fastp json
    with open(fastp_json, 'r') as f:
        data = json.loads(f.read())
        n_reads_tot = int(data['summary']['before_filtering']['total_reads'] / 2)

    # grab aligned read numbers from bam
    n_reads_mRNA = get_read_number(mRNA_bam)
    n_reads_mRNA_ded = get_read_number(mRNA_bam_ded)
    n_reads_rRNA = get_read_number(rRNA_bam)
    n_reads_rRNA_ded = get_read_number(rRNA_bam_ded)

    ## DMS stats------------------------------------------------------------------------------------
    
    # get signal-to-noise
    bin_size = 0.05
    genome = Genome(genome_fasta)
    genome_rrna = Genome(genome_fasta, coords=(11, 451575, 458433), reverse=True)
    dms_mRNA = DMS_Profile(dms_mRNA_pkl, 'blank', genome, with_stats=True)
    dms_rRNA = Targeted_DMS_Profile(dms_rRNA_pkl, 'blank', genome_rrna, reverse=True, with_stats=True)

    s2n_mRNA, x_mRNA, y_gu_mRNA, y_ac_mRNA = get_s2n(dms_mRNA, genome, bin_size=bin_size)
    s2n_rRNA, x_rRNA, y_gu_rRNA, y_ac_rRNA = get_s2n(dms_rRNA, genome_rrna, bin_size=bin_size, targeted=True)

    # plotting
    ax2.bar(x_mRNA, y_ac_mRNA, bin_size, color='mediumseagreen')
    ax2.step(x_mRNA, y_gu_mRNA, bin_size, color='r')
    ax2.set_title(f"mRNA, s2n = {s2n_mRNA:.2f}, n_AC = {np.sum(np.hstack(dms_mRNA.mask_p) & np.hstack(genome.mask_seq)) + np.sum(np.hstack(dms_mRNA.mask_m) & ~np.hstack(genome.mask_seq))}")
    ax2.set_xlabel('log10(mutation rate)')

    ax3.bar(x_rRNA, y_ac_rRNA, bin_size, color='mediumseagreen', label='A/C')
    ax3.step(x_rRNA, y_gu_rRNA, bin_size, color='r', label='G/U')
    ax3.set_title(f"rRNA, s2n = {s2n_rRNA:.2f}, n_AC = {np.sum(dms_rRNA.mask & genome_rrna.mask_seq)}")
    ax3.set_xlabel('log10(mutation rate)')
    ax3.legend()


    ## Insert size distribution-------------------------------------------------------------------

    # get data
    len_mRNA = PRO_Profile(pro_mRNA_pkl, 'blank', genome, with_stats=True).stats

    # plotting
    ax1.hist(np.array(len_mRNA), bins=np.arange(0, 600, 1), density=False, alpha=0.6)
    ax1.set_xlabel('insert size')
    ax1.set_title(f"median length: {np.median(len_mRNA):.2f}")


    ## DMS read stats----------------------------------------------------------------------------
    n_reads_mRNA_dms = np.sum([len(i) for i in dms_mRNA.stats])
    n_reads_rRNA_dms = np.sum([len(i) for i in dms_rRNA.stats])

    ax4.hist([item for sublist in dms_mRNA.stats for item in sublist], bins=np.arange(0, 15, 1), align='left', density=True)
    ax4.set_xticks(np.arange(0, 15, 1))
    ax4.set_xlabel('mutations per read')
    ax4.set_title('mRNA')

    ax5.hist([item for sublist in dms_rRNA.stats for item in sublist], bins=np.arange(0, 15, 1), align='left', density=True)
    ax5.set_xticks(np.arange(0, 15, 1))
    ax5.set_xlabel('mutations per read')
    ax5.set_title('rRNA')

    ## example profiles----------------------------------------------------------------------
    chrom = 6
    coords = (882812, 882812+200)

    dms_mRNA.plot_profile(chrom, '-', coords[0], coords[1], genome, ax6, cmap_loc='/gpfs/gibbs/project/neugebauer/ls2286/lsStruct/scripts/cmap.txt')
    ax6.set_title('mRNA: TDH3')
    dms_rRNA.plot_profile(3252, 3252+200, genome_rrna, ax7, cmap_loc='/gpfs/gibbs/project/neugebauer/ls2286/lsStruct/scripts/cmap.txt')
    plt.setp([ax6, ax7], ylim=(0, 0.15))
    ax7.set_title('rRNA')

    plt.savefig(output_path + '/' + sample_name + '_QC.pdf')

    ## save text log-----------------------------------------------------------------------
    log_df = pd.DataFrame.from_dict({'raw': n_reads_tot,
                        'mRNA': n_reads_mRNA,
                        'rRNA': n_reads_rRNA,
                        'mRNA deduplicated': n_reads_mRNA_ded, 
                        'rRNA deduplicated': n_reads_rRNA_ded,
                        'mRNA passed DMS filters': n_reads_mRNA_dms,
                        'rRNA passed DMS filters': n_reads_rRNA_dms}, orient='index').reset_index()

    log_df['ratios'] = [100, 100*n_reads_mRNA/n_reads_tot, 100*n_reads_rRNA/n_reads_tot, 100*n_reads_mRNA_ded/n_reads_mRNA, 100*n_reads_rRNA_ded/n_reads_rRNA, 100*n_reads_mRNA_dms/n_reads_mRNA_ded, 100*n_reads_rRNA_dms/n_reads_rRNA_ded]

    log_df.to_csv(output_path + '/' + sample_name + '.log', header=None, index=None)

def get_s2n(dms, genome, bin_size=0.05, targeted=False):

    if targeted:
        signa = dms.rat[dms.mask & genome.mask_seq]
        noise = dms.rat[dms.mask & ~genome.mask_seq]
    else:
        signa = np.hstack([np.hstack(dms.rat_p)[np.hstack(dms.mask_p) & np.hstack(genome.mask_seq)], np.hstack(dms.rat_m)[np.hstack(dms.mask_m) & ~np.hstack(genome.mask_seq)]])
        noise = np.hstack([np.hstack(dms.rat_p)[np.hstack(dms.mask_p) & ~np.hstack(genome.mask_seq)], np.hstack(dms.rat_m)[np.hstack(dms.mask_m) & np.hstack(genome.mask_seq)]])

    # get the distributions
    x = np.arange(-4, 0, bin_size)[:-1] + 0.5*bin_size
    y_gu = np.histogram(np.log10(noise), bins=np.arange(-4, 0, bin_size), density=True)[0]
    y_ac = np.histogram(np.log10(signa), bins=np.arange(-4, 0, bin_size), density=True)[0]
    
    # calculate s2n
    snr = np.median(signa*dms.fac_norm) / np.median(noise*dms.fac_norm)
    
    return(snr, x, y_gu, y_ac)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('sample_name')
    parser.add_argument('output_path')
    parser.add_argument('fastp_json')
    parser.add_argument('mRNA_bam')
    parser.add_argument('mRNA_bam_ded')
    parser.add_argument('rRNA_bam')
    parser.add_argument('rRNA_bam_ded')
    parser.add_argument('dms_mRNA_pkl')
    parser.add_argument('dms_rRNA_pkl')
    parser.add_argument('genome_fasta')
    parser.add_argument('pro_mRNA_pkl')

    args = parser.parse_args()

    QC(args.sample_name, args.output_path, args.fastp_json, args.mRNA_bam, args.mRNA_bam_ded, args.rRNA_bam, args.rRNA_bam_ded, args.dms_mRNA_pkl, args.dms_rRNA_pkl, args.genome_fasta, args.pro_mRNA_pkl)
