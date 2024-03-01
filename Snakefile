## load locations of tools from config file (just for readability of bash code)
BWA=config['BWA']
SAMTOOLS=config['SAMTOOLS']
DEDUPBAM=config['DEDUPBAM']
PYTHON=config['PYTHON']
FASTP=config['FASTP']
STAR=config['STAR']
UMIDEDUP=config['UMIDEDUP']
DESEQ=config['DESEQ']
SPLICEQ=config['SPLICEQ']
HDP=config['HDP']

## collect sample information and gene lists for targeted analysis
import pandas as pd
import numpy as np

samples = pd.read_csv(config["SAMPLES"], delimiter=',').set_index('sample_name')
cotrx = pd.read_csv(config["COTRX_FILE"], delimiter=',').set_index('name')
read2_regions = pd.read_csv(config["READ2_COORDS"], delimiter=',').set_index('name')

def get_output(samples, OUTPUTPATH, RUN_DMS_AGG, RUN_DMS_R2, RUN_R2_REGIONS, RUN_DMS_COTRX_RRNA, RUN_DMS_COTRX, RUN_PRO_SIGNAL, RUN_PRO_SIGNAL_RRNA, RUN_DESEQ, RUN_HDPROBE, DESEQ_CTRL, RUN_SPLICEQ, RUN_QC):
    # preprocessing
    if RUN_QC:
        out_files = [OUTPUTPATH + i + '.log' for i in list(samples.index)]
    else:
        out_files = []
    # DMS_aggregate
    if RUN_DMS_AGG:
        for i in [OUTPUTPATH + 'dms/' + i + '_mRNA_agg.pkl' for i in list(samples.index)]:
            out_files.append(i)
        for i in [OUTPUTPATH + 'dms/' + i + '_rRNA_agg.pkl' for i in list(samples.index)]:
            out_files.append(i)
    # DMS_read2
    if RUN_DMS_R2:
        for i in [OUTPUTPATH + 'dms/' + i + '_mRNA_read2.pkl' for i in list(samples.index)]:
            out_files.append(i)
        for i in [OUTPUTPATH + 'dms/' + i + '_rRNA_read2.pkl' for i in list(samples.index)]:
            out_files.append(i)
    # DMS_cotrx_rRNA
    if RUN_DMS_COTRX_RRNA:
        for i in [OUTPUTPATH + 'dms/' + i + '_rRNA.pkl' for i in list(samples.index)]:
            out_files.append(i)
    # PRO_signal_mRNA
    if RUN_PRO_SIGNAL:
        for i in [OUTPUTPATH + 'pro/' + i + '_mRNA_pro.pkl' for i in list(samples.index)]:
            out_files.append(i)
    # PRO_signal_rRNA
    if RUN_PRO_SIGNAL_RRNA:
        for i in [OUTPUTPATH + 'pro/' + i + '_rRNA_pro.pkl' for i in list(samples.index)]:
            out_files.append(i)
    # cotranscriptional targets that are not rRNA
    if RUN_DMS_COTRX:
        for j in list(cotrx.index):
            for i in [OUTPUTPATH + 'dms/' + i + '_cotrx_' + j + '.pkl' for i in list(samples.index)]:
                out_files.append(i)
    # read2 meta for specific regions
    if RUN_R2_REGIONS:
        for j in list(read2_regions.index):
            for i in [OUTPUTPATH + 'dms/' + i + '_read2_' + j + '.pkl' for i in list(samples.index)]:
                out_files.append(i)
    # differential expression
    if RUN_DESEQ:
        out_files.append(OUTPUTPATH + 'deseq/deseq_results.pkl')
    # HDProbe
    if RUN_HDPROBE:
        for i in list(np.unique([j.split('_')[1] for j in samples.index if j.split('_')[1] != DESEQ_CTRL])):
            out_files.append(OUTPUTPATH + 'HDProbe/' + i + '_mrna.csv')
            out_files.append(OUTPUTPATH + 'HDProbe/' + i + '_rrna.csv')
    # SPLICE-q
    if RUN_SPLICEQ:
        for i in [OUTPUTPATH + 'splice/' + i + '_SE.csv' for i in list(samples.index)]:
            out_files.append(i)
    return(out_files)   



## run rule
rule all:
    input:
        get_output(samples, config['OUTPUTPATH'], config['RUN_DMS_AGG'], config['RUN_DMS_R2'], config['RUN_R2_REGIONS'], config['RUN_DMS_COTRX_RRNA'], config['RUN_DMS_COTRX'], config['RUN_PRO_SIGNAL'], config['RUN_PRO_SIGNAL_RRNA'], config['RUN_DESEQ'], config['RUN_HDPROBE'], config['DESEQ_CTRL'], config['RUN_SPLICEQ'], config['RUN_QC'])


## modules
rule fastp:
    input:
        r1=ancient(config['DATAPATH'] + "{sample}_R1.fastq.gz"),
        r2=ancient(config['DATAPATH'] + "{sample}_R2.fastq.gz")
    output:
        r1=config['OUTPUTPATH'] + "fastq/{sample}_R1_fp.fastq.gz",
        r2=config['OUTPUTPATH'] + "fastq/{sample}_R2_fp.fastq.gz",
        html=config['OUTPUTPATH'] + "fastq/{sample}_fastp.html",
        json=config['OUTPUTPATH'] + "fastq/{sample}_fastp.json"
    threads: config['THREADS']
    params:
        l=config['MINREADLENGTH'],
        umi_len=config['UMILENGTH'],
        umi_prefix=config['UMIPREFIX'],
        trimfront=config['TRIMFRONT2'],
        adapter_r1 = lambda wildcards: samples.loc[wildcards.sample, 'adapter_r1'],
        adapter_r2 = lambda wildcards: samples.loc[wildcards.sample, 'adapter_r2']
    shell:
        "{FASTP} -l {params.l} -c --trim_poly_g --adapter_sequence {params.adapter_r1} --adapter_sequence_r2 {params.adapter_r2} "
        "--umi --umi_loc read1 --umi_len {params.umi_len} --umi_prefix {params.umi_prefix} --dedup "
        "--trim_front2 {params.trimfront} -w {threads} --json {output.json} --html {output.html} "
        "-i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2}"

rule bwa_index:
    input:
        ancient(config['GENOMEFASTA'])
    output:
        config['GENOMEFASTA'] + ".sa"
    shell:
        "{BWA} index {input}"


rule bwa_mem:
    input:
        sample=[config['OUTPUTPATH'] + "fastq/{sample}_R1_fp.fastq.gz", config['OUTPUTPATH'] + "fastq/{sample}_R2_fp.fastq.gz"],
        requir=config['GENOMEFASTA'] + ".sa",
        genome=config['GENOMEFASTA']
    output:
        temp(config['OUTPUTPATH'] + "bam/{sample}.bam")
    params:
        m="500M",
        L="5,999" if config['CLIPPENALTY_5'] else "5,5"
    threads: config['THREADS']
    shell:
        "{BWA} mem -t {threads} -L {params.L} {input.genome} {input.sample} | "
        "{SAMTOOLS} sort --threads {threads} -m {params.m} -n -o {output}"

rule star_index:
    input:
        fasta=ancient(config["GENOMEFASTA"]), gtf=ancient(config["GENOMEANN"])
    output:
        directory(config["STARGENOMEDIR"])
    params:
        genomedir=config["STARGENOMEDIR"]
    threads: 2
    shell:
        # for yeast
        "{STAR} --runThreadN 1 --runMode genomeGenerate --genomeDir {params.genomedir} --genomeFastaFiles {input.fasta} --sjdbGTFfile {input.gtf} --genomeSAindexNbases 10 --sjdbGTFtagExonParentTranscript Parent --sjdbGTFfeatureExon CDS"
        # for mouse
        #"{STAR} --runThreadN {threads} --runMode genomeGenerate --genomeDir {params.genomedir} --genomeFastaFiles {input.fasta} --sjdbGTFfile {input.gtf}"


rule STAR:
    input:
        samples=[config["OUTPUTPATH"] + "fastq/{sample}_R1_fp.fastq.gz", config["OUTPUTPATH"] + "fastq/{sample}_R2_fp.fastq.gz"],
        genome=config["STARGENOMEDIR"]
    output:
        starout=temp(config["OUTPUTPATH"] + "bam/{sample}Aligned.out.bam"),
        samout=temp(config["OUTPUTPATH"] + "bam/{sample}_star.bam")
    threads: config["THREADS"]
    params:
        l=config['MINREADLENGTH'],
        m="500M",
        odir=config["OUTPUTPATH"],
        max_insert=config["STAR_MAX_INSERT"],
        intron_min=config["STAR_INTRON_MIN"],
        intron_max=config["STAR_INTRON_MAX"],
        splice_OH=config["STAR_SJD_OHANG"]
    shell:
        """
        {STAR} --runThreadN {threads} --alignSJoverhangMin {params.splice_OH} --alignMatesGapMax {params.max_insert} --alignIntronMin {params.intron_min} --alignIntronMax {params.intron_max} --readFilesCommand zcat --outSAMtype BAM Unsorted --outFileNamePrefix {params.odir}bam/{wildcards.sample} --genomeDir {input.genome} --readFilesIn {input.samples}
        {SAMTOOLS} sort --threads {threads} -m {params.m} -n -o {output.samout} {output.starout}
        """


rule samtools_fixmate:
    input:
        config['OUTPUTPATH'] + "bam/{sample}_star.bam" if config["USE_STAR"] else config['OUTPUTPATH'] + "bam/{sample}.bam"
    output:
        temp(config['OUTPUTPATH'] + "bam/{sample}_fixed.bam")
    threads: config['THREADS']
    shell:
        "{SAMTOOLS} fixmate --threads {threads} -r -m {input} {output}"


rule samtools_sort:
    input:
        bam=config['OUTPUTPATH'] + "bam/{sample}_fixed.bam",
        genome=config['GENOMEFASTA']
    output:
        bam=temp(config['OUTPUTPATH'] + "bam/{sample}_fixed_sorted.bam"),
        bai=temp(config['OUTPUTPATH'] + "bam/{sample}_fixed_sorted.bam.bai")
    params:
        m="500M",
        star=config['USE_STAR']
    threads: config['THREADS']
    shell:
        """
        if [ {params.star} = True ]; then
            echo "Calling MD tag for STAR alignments."
            {SAMTOOLS} sort --threads {threads} -m {params.m} {input.bam} | {SAMTOOLS} calmd -b - {input.genome} > {output.bam}
        else
            echo "Sorting, MD tag already there from bwa."
            {SAMTOOLS} sort --threads {threads} -m {params.m} -o {output.bam} {input.bam}
        fi
        {SAMTOOLS} index {output.bam}
        """

rule separate_rRNA:
    input:
        sample=config['OUTPUTPATH'] + "bam/{sample}_fixed_sorted.bam",
        rcoord=config['RRNABED']
    output:
        rRNA_bam=config['OUTPUTPATH'] + "bam/{sample}_rRNA.bam",
        rRNA_bai=config['OUTPUTPATH'] + "bam/{sample}_rRNA.bam.bai",
        mRNA_bam=temp(config['OUTPUTPATH'] + "bam/{sample}_fin.bam")
    shell:
        """
        {SAMTOOLS} view -o {output.rRNA_bam} -U {output.mRNA_bam} -L {input.rcoord} {input.sample}
        {SAMTOOLS} index {output.rRNA_bam}
        """

rule filter_mRNA:
    input:
        config['OUTPUTPATH'] + "bam/{sample}_fin.bam"
    output:
        bam=config['OUTPUTPATH'] + "bam/{sample}_fin_MQ.bam",
        bai=config['OUTPUTPATH'] + "bam/{sample}_fin_MQ.bam.bai"
    params:
        q=config['MQFILTER']
    shell:
        """
        {SAMTOOLS} view -q {params.q} -o {output.bam} {input}
        {SAMTOOLS} index {output.bam}
        """

rule deduplicate_mRNA:
    input:
        config['OUTPUTPATH'] + "bam/{sample}_fin_MQ.bam"
    output:
        bam=config['OUTPUTPATH'] + "bam/{sample}_fin_dedup.bam",
        bai=config['OUTPUTPATH'] + "bam/{sample}_fin_dedup.bam.bai"
    params:
        umi=config['UMIPREFIX'] + '_',
        umi_len=config['UMILENGTH'],
        umi_bool= ' --without_umi ' if config['WITHOUTUMI'] else ' ',
        seq_bool= ' --ignore_seq ' if config['IGNORESEQ'] else ' '
    resources:
        mem_mb=config['DEDUPMEM']
    shell:
        """
        #{PYTHON} {DEDUPBAM} -umi {params.umi}{params.seq_bool}-umi_len {params.umi_len}{params.umi_bool}-o {output.bam} {input}
        {UMIDEDUP} bam -i {input} -o {output.bam} --umi-sep _ --paired --two-pass --merge avgqual --remove-unpaired --remove-chimeric
        {SAMTOOLS} index {output.bam}
        """

rule deduplicate_rRNA:
    input:
        config['OUTPUTPATH'] + "bam/{sample}_rRNA.bam"
    output:
        bam=config['OUTPUTPATH'] + "bam/{sample}_rRNA_dedup.bam",
        bai=config['OUTPUTPATH'] + "bam/{sample}_rRNA_dedup.bam.bai"
    params:
        umi=config['UMIPREFIX'] + '_',
        umi_len=config['UMILENGTH'],
        umi_bool= ' --without_umi ' if config['WITHOUTUMI'] else ' ',
        seq_bool= ' '
    resources:
        mem_mb=config['DEDUPMEM']
    shell:
        """
        #{PYTHON} {DEDUPBAM} -umi {params.umi}{params.seq_bool}-umi_len {params.umi_len}{params.umi_bool}-o {output.bam} {input}
        {UMIDEDUP} bam -i {input} -o {output.bam} --umi-sep _ --paired --two-pass --merge avgqual --remove-unpaired --remove-chimeric
        {SAMTOOLS} index {output.bam}
        """

rule QC:
    input:
        path_out  = config['OUTPUTPATH'],
        fastp_json = config['OUTPUTPATH'] + "fastq/{sample}_fastp.json",
        mRNA_bam   = config['OUTPUTPATH'] + "bam/{sample}_fin_MQ.bam",
        mRNA_bam_ded = config['OUTPUTPATH'] + "bam/{sample}_fin_dedup.bam",
        rRNA_bam = config['OUTPUTPATH'] + "bam/{sample}_rRNA.bam",
        rRNA_bam_ded = config['OUTPUTPATH'] + "bam/{sample}_rRNA_dedup.bam",
        dms_mRNA_pkl = config['OUTPUTPATH'] + "dms/{sample}_mRNA_agg.pkl",
        dms_rRNA_pkl = config['OUTPUTPATH'] + "dms/{sample}_rRNA_agg.pkl",
        genome_fasta = config["GENOMEFASTA"],
        pro_mRNA_pkl = config['OUTPUTPATH'] + "pro/{sample}_mRNA_pro.pkl"
    output:
        config['OUTPUTPATH'] + "{sample}.log",
        config['OUTPUTPATH'] + "{sample}_QC.pdf"
    shell:
        "{PYTHON} ../../scripts/run_QC.py {wildcards.sample} {input.path_out} {input.fastp_json} {input.mRNA_bam} {input.mRNA_bam_ded} {input.rRNA_bam} {input.rRNA_bam_ded} {input.dms_mRNA_pkl} {input.dms_rRNA_pkl} {input.genome_fasta} {input.pro_mRNA_pkl}"

rule DMS_aggregate:
    input:
        bam=config['OUTPUTPATH'] + "bam/{sample}_fin_dedup.bam" if config['DEDUP_MRNA'] else config['OUTPUTPATH'] + "bam/{sample}_fin_MQ.bam",
        genome=config['GENOMEFASTA']
    output:
        config['OUTPUTPATH'] + "dms/{sample}_mRNA_agg.pkl"
    threads: config['THREADS_DMS']
    params:
        clip_filter=config['CLIPFILTER'],
        clip_3=config['CLIP_3'],
        clip_filter_3=config['FILTER_3_CLIPPED'],
        allow_polyA=config['ALLOW_POLYA']
    run:
        from src.dms_utils import coverage_mutations_PE
        import pickle

        mut_p, cov_p, mut_m, cov_m, stats = coverage_mutations_PE(input.bam, input.genome, return_read_stats=True,  n_cpu=threads, clip_filter=params.clip_filter, clip_3=params.clip_3, clip_filter_3=params.clip_filter_3, allow_polyA_clipping=params.allow_polyA)
        with open(output[0], 'wb') as f:
            pickle.dump((mut_p, cov_p, mut_m, cov_m, stats), f)


rule DMS_aggregate_rRNA:
    input:
        bam=config['OUTPUTPATH'] + "bam/{sample}_rRNA_dedup.bam"  if config['DEDUP_RRNA'] else config['OUTPUTPATH'] + "bam/{sample}_rRNA.bam",
        genome=config['GENOMEFASTA']
    output:
        config['OUTPUTPATH'] + "dms/{sample}_rRNA_agg.pkl"
    threads: 1
    params:
        clip_filter=config['CLIPFILTER'],
        clip_3=config['CLIP_3']
    run:
        from src.dms_utils import coverage_mutations_PE
        import pickle

        _, _, mut_m, cov_m, stats = coverage_mutations_PE(input.bam, input.genome, return_read_stats=True, n_cpu=threads, clip_filter=params.clip_filter, clip_3=params.clip_3, clip_filter_3=True)
        with open(output[0], 'wb') as f:
            pickle.dump((mut_m[11][451575:458433], cov_m[11][451575:458433], stats), f)


rule DMS_read2_mRNA:
    input:
        bam=config['OUTPUTPATH'] + "bam/{sample}_fin_dedup.bam" if config['DEDUP_MRNA'] else config['OUTPUTPATH'] + "bam/{sample}_fin_MQ.bam",
        genome=config['GENOMEFASTA']
    output:
        config['OUTPUTPATH'] + "dms/{sample}_mRNA_read2.pkl"
    params:
        length=config['READ2DMSLEN']
    threads: config['THREADS_DMS']
    run:
        from src.pro_utils import DMS_read2_parallel
        import pickle

        r2_mut_dms, r2_cov_dms = DMS_read2_parallel(input.bam, params.length, input.genome, n_cpu=threads)
        with open(output[0], 'wb') as f:
            pickle.dump((r2_mut_dms, r2_cov_dms), f)


rule DMS_read2_rRNA:
    input:
        bam=config['OUTPUTPATH'] + "bam/{sample}_rRNA_dedup.bam"  if config['DEDUP_RRNA'] else config['OUTPUTPATH'] + "bam/{sample}_rRNA.bam",
        genome=config['GENOMEFASTA']
    output:
        config['OUTPUTPATH'] + "dms/{sample}_rRNA_read2.pkl"
    params:
        length=config['READ2DMSLEN']
    run:
        from src.pro_utils import DMS_read2_parallel
        import pickle

        # not running parallel because rRNA isn't currently parallelized (single chromosome) sue me
        r2_mut_dms, r2_cov_dms = DMS_read2_parallel(input.bam, params.length, input.genome, n_cpu=1)
        with open(output[0], 'wb') as f:
            pickle.dump((r2_mut_dms, r2_cov_dms), f)

rule DMS_read2_region:
    input:
        bam=config['OUTPUTPATH'] + "bam/{sample}_fin_dedup.bam" if config['DEDUP_MRNA'] else config['OUTPUTPATH'] + "bam/{sample}_fin_MQ.bam",
        genome=config['GENOMEFASTA'],
        regions_file=lambda wildcards: read2_regions.loc[f"{wildcards.region}", "location"]
    output:
        config['OUTPUTPATH'] + "dms/{sample}_read2_{region}.pkl"
    params:
        length=config['READ2DMSLEN'],
        clipfilter=config['CLIPFILTER']
    threads: config['THREADS_DMS']
    run:
        from src.pro_utils import DMS_read2_region_parallel
        import pickle

        r2_mut_dms, r2_cov_dms = DMS_read2_region_parallel(input.bam, params.length, input.genome, input.regions_file, clip_filter=params.clipfilter, n_cpu=threads)
        with open(output[0], 'wb') as f:
            pickle.dump((r2_mut_dms, r2_cov_dms), f)


rule DMS_cotrx_rRNA:
# this is currently hard-coded for yeast!
    input:
        config['OUTPUTPATH'] + "bam/{sample}_rRNA_dedup.bam"  if config['DEDUP_RRNA'] else config['OUTPUTPATH'] + "bam/{sample}_rRNA.bam"
    output:
        config['OUTPUTPATH'] + "dms/{sample}_rRNA.pkl"
    params:
        clipfilter=config['CLIPFILTER']
    run:
        from src.pro_utils import PRO_DMS_signal
        from src.CotrxMatrix import sparse_to_dict
        import pickle
        import numpy as np

        mat_mut, mat_cov = PRO_DMS_signal(input[0], 'chrXII', 451575, 458432+1, '-', clip_filter=params.clipfilter)

        # convert to dict to reduce disk space when saving
        mut_d = sparse_to_dict(mat_mut)
        cov_d = sparse_to_dict(mat_cov)

        with open(output[0], 'wb') as f:
            pickle.dump((mut_d, cov_d, mat_mut.shape), f)

rule DMS_cotrx:
    input:
        config['OUTPUTPATH'] + "bam/{sample}_fin_dedup.bam"  if config['DEDUP_MRNA'] else config['OUTPUTPATH'] + "bam/{sample}_fin_MQ.bam"
    output:
        config['OUTPUTPATH'] + "dms/{sample}_cotrx_{cotrxname}.pkl"
    params:
        ref = lambda wildcards: cotrx.loc[wildcards.cotrxname, 'ref'],
        strand = lambda wildcards: cotrx.loc[wildcards.cotrxname, 'strand'],
        start = lambda wildcards: cotrx.loc[wildcards.cotrxname, 'start'],
        end = lambda wildcards: cotrx.loc[wildcards.cotrxname, 'end'],
        clipfilter=config['CLIPFILTER']
    run:
        from src.pro_utils import PRO_DMS_signal
        from src.CotrxMatrix import sparse_to_dict
        import pickle

        mat_mut, mat_cov = PRO_DMS_signal(input[0], params.ref, params.start, params.end, params.strand, clip_filter=params.clipfilter)

        # convert to dict to reduce disk space when saving
        mut_d = sparse_to_dict(mat_mut)
        cov_d = sparse_to_dict(mat_cov)

        with open(output[0], 'wb') as f:
            pickle.dump((mut_d, cov_d, mat_mut.shape), f)


rule PRO_signal_mRNA:
    input:
        bam=config['OUTPUTPATH'] + "bam/{sample}_fin_dedup.bam" if config['DEDUP_MRNA'] else config['OUTPUTPATH'] + "bam/{sample}_fin_MQ.bam",
        genome=config['GENOMEFASTA']
    output:
        config['OUTPUTPATH'] + "pro/{sample}_mRNA_pro.pkl"
    run:
        from src.pro_utils import map_read_ends
        import pickle

        polym_p, _, polym_m, _, l = map_read_ends(input.bam, input.genome, with_lengths=True)

        with open(output[0], 'wb') as f:
            pickle.dump((polym_p, polym_m, l), f)

rule PRO_signal_rRNA:
    input:
        bam=config['OUTPUTPATH'] + "bam/{sample}_rRNA_dedup.bam"  if config['DEDUP_RRNA'] else config['OUTPUTPATH'] + "bam/{sample}_rRNA.bam",
        genome=config['GENOMEFASTA']
    output:
        config['OUTPUTPATH'] + "pro/{sample}_rRNA_pro.pkl"
    run:
        from src.pro_utils import map_read_ends
        import pickle

        polym_p, _, polym_m, _, l = map_read_ends(input.bam, input.genome, with_lengths=True)

        with open(output[0], 'wb') as f:
            pickle.dump((polym_p, polym_m, l), f)


rule featureCounts:
    input:
        bams=expand(config['OUTPUTPATH'] + "bam/{sample}_fin_dedup.bam", sample=list(samples.index)) if config['DEDUP_MRNA'] else expand(config['OUTPUTPATH'] + "bam/{sample}_fin_MQ.bam", sample=list(samples.index)),
        ann=config['GENOMEANN_gtf']
    output:
        config['OUTPUTPATH'] + "deseq/diff.counts.tsv"
    threads: config['THREADS_DMS']
    shell:
        # yeast
        "featureCounts -T {threads} -Q 10 --primary -s 1 --countReadPairs -pBP -t CDS -g transcript_id -a {input.ann} -o {output} {input.bams}"
        # mouse
        #"featureCounts -T {threads} -Q 10 --primary -s 1 -pC -t exon -g gene_id --extraAttributes gene_name -a {input.ann} -o {output} {input.bams}"

rule DEseq2:
    input:
        config['OUTPUTPATH'] + "deseq/diff.counts.tsv"
    output:
        config['OUTPUTPATH'] + "deseq/deseq_results.pkl"
    params:
        ctrl=config['DESEQ_CTRL']
    threads: config['THREADS_DMS']
    shell:
        "{PYTHON} {DESEQ} -ctrl {params.ctrl} -c {threads} -o {output} {input}"


rule HDProbe_prep:
    input:
        dms_data=expand(config['OUTPUTPATH'] + "dms/{sample}_mRNA_agg.pkl", sample=list(samples.index)),
        fasta=config["GENOMEFASTA"]
    params:
        sample_names=list(samples.index),
        cutoff=config['MIN_COV']
    output:
        hdp_in=config['OUTPUTPATH'] + "HDProbe/mut_cov_mrna.tsv"
    run:
        from src.DMS_Profile import write_HDProbe, DMS_Profile, Genome

        genome = Genome(input.fasta)
        write_HDProbe(output.hdp_in, [DMS_Profile(i, name, genome) for i, name in zip(input.dms_data, params.sample_names)], genome.seq, min_cov=params.cutoff)

rule HDProbe:
    input:
        config['OUTPUTPATH'] + "HDProbe/mut_cov_mrna.tsv"
    output:
        config['OUTPUTPATH'] + "HDProbe/{condition}_mrna.csv"
    params:
        control=config['DESEQ_CTRL'],
        cutoff=700
    resources:
        mem_mb=16000
    threads: 1
    priority: -50
    shell:
        """
        conda deactivate
    module purge
        module load R-bundle-Bioconductor/3.15-foss-2020b-R-4.2.0
        Rscript ../../scripts/run_HDP.R -a {params.control} -b {wildcards.condition} --hdploc {HDP} --mutcov {input} --outpath {output} --cutoff {params.cutoff}
        """

rule HDProbe_prep_rrna:
    input:
        dms_data=expand(config['OUTPUTPATH'] + "dms/{sample}_rRNA_agg.pkl", sample=list(samples.index)),
        fasta=config["GENOMEFASTA"]
    params:
        sample_names=list(samples.index),
        cutoff=config['MIN_COV']
    output:
        hdp_in=config['OUTPUTPATH'] + "HDProbe/mut_cov_rrna.tsv"
    run:
        from src.DMS_Profile import targeted_write_HDProbe, Targeted_DMS_Profile, Genome

        genome = Genome(input.fasta, coords=(11, 451575, 458433), reverse=True)
        targeted_write_HDProbe(output.hdp_in, [Targeted_DMS_Profile(i, name, genome, reverse=True) for i, name in zip(input.dms_data, params.sample_names)], genome.seq, min_cov=params.cutoff)

rule HDProbe_rrna:
    input:
        config['OUTPUTPATH'] + "HDProbe/mut_cov_rrna.tsv"
    output:
        config['OUTPUTPATH'] + "HDProbe/{condition}_rrna.csv"
    params:
        control=config['DESEQ_CTRL'],
        cutoff=700
    threads: 1
    resources:
        mem_mb=16000
    priority: -50
    shell:
        """
        module purge
        module load R-bundle-Bioconductor/3.15-foss-2020b-R-4.2.0
        Rscript ../../scripts/run_HDP.R -a {params.control} -b {wildcards.condition} --hdploc {HDP} --mutcov {input} --outpath {output} --cutoff {params.cutoff}
        """

rule Spliceq:
    input:
        bam=config['OUTPUTPATH'] + "bam/{sample}_fin_dedup.bam" if config['DEDUP_MRNA'] else config['OUTPUTPATH'] + "bam/{sample}_fin_MQ.bam",
        genome=config['GENOMEANN_gtf']
    output:
        config['OUTPUTPATH'] + "splice/{sample}_SE.csv"
    shell:
        "{PYTHON} {SPLICEQ} -c 10 -b {input.bam} -g {input.genome} -o {output} --quiet"











