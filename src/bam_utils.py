from collections import defaultdict
import re
import subprocess

def read_pair_generator(bam, reference=None, start=None, end=None):
    """
    Generate read pairs in a BAM file or within a region string.
    from https://www.biostars.org/p/306041/
    
    Args:
        bam (pysam.AlignmentFile): BAM file.
        reference (str): Reference sequence name.
        start (int): Start position for region filtering.
        end (int): End position for region filtering.

    Yields:
        Tuple[pysam.AlignedSegment, pysam.AlignedSegment]: Paired reads.
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(reference=reference, start=start, end=end):
        if not read.is_proper_pair or read.is_secondary or read.is_supplementary or read.is_unmapped:
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]


def read_ends(read):
    """
    Get the mapped ends of a read.

    Args:
        read (pysam.AlignedSegment): Input read.

    Returns:
        Tuple[str, int, int]: Reference name, mapped 5' end, mapped 3' end.
    """
    if read.is_unmapped or read.is_secondary or read.is_supplementary or not read.is_proper_pair:
        return(None)

    if read.is_reverse:
        mapped_5 = read.reference_end-1
        mapped_3 = read.reference_start
    else:
        mapped_5 = read.reference_start
        mapped_3 = read.reference_end-1


    return(read.reference_name, mapped_5, mapped_3)


def get_clipped(read):
    """
    Get the clipped bases from a read.

    Args:
        read (pysam.AlignedSegment): Input read.

    Returns:
        Tuple[str, str]: Clipped 3' end and clipped 5' end.
    """

    if read.is_reverse:
        reverse = True
        read_seq = reverse_complement(read.get_forward_sequence()).upper()
    else:
        reverse = False
        read_seq = read.get_forward_sequence().upper()  # Sequence of the read

    clipped_3 = ''
    clipped_5 = ''

    i = 0  # Pos in the ref sequence
    j = 0  # Pos in the read sequence
    l = 0  # Pos in the ref position list
    cigar =re.findall(r'(\d+)([A-Z]{1})', read.cigarstring)
    op_index = 0
    while op_index < len(cigar):  # Each CIGAR operation
        op = cigar[op_index]
        desc, length = op[1], int(op[0])

        if desc == 'M':  # Match or mismatch
            i += length  # Update ref index
            j += length  # Update read index
            l += length  # Update position index

        elif desc == 'D': # Deletion
            i += length  # Update ref index

        elif desc == 'I':  # Insertion
            j += length  # Update read index
            l += length  

        elif desc == 'S':  # Soft clipping
            if (not reverse and op_index == len(cigar) - 1) or (reverse and op_index == 0):  # Soft clipped at the 3'-end
                l += length  # Update position index (soft-clipped bases are part of the position list)
                clipped_3 = read_seq[j:j+length] # clipped 3'-end bases for studying poly(A) tail properties
            else: # Soft clipped at 5'-end
                clipped_5 = read_seq[j:j+length] # clipped 5'-end bases for studying RT-stop properties
                l += length
            j += length  # Update read index

        op_index += 1

    return(clipped_3, clipped_5)


def complement(dna):
    """
    Return the complement of a DNA sequence.

    Args:
        dna (str): DNA sequence.

    Returns:
        str: Complemented DNA sequence.
    """

    def _complement(x):
        if x == 'A':
            return('T')
        elif x == 'T':
            return('A')
        elif x == 'C':
            return('G')
        elif x == 'G':
            return('C')
        elif x == 'U':
            return('A')
        elif x == 'a':
            return('t')
        elif x == 't':
            return('a')
        elif x == 'c':
            return('g')
        elif x == 'g':
            return('c')
        elif x == 'u':
            return('a')
        elif x == 'N':
            return('N')
    return(''.join([_complement(x) for x in dna]))

def reverse_complement(dna):
    """
    Return the reverse complement of a DNA sequence.

    Args:
        dna (str): DNA sequence.

    Returns:
        str: Reverse complemented DNA sequence.
    """

    return(complement(dna)[::-1])

def get_read_number(bamfile, samtools_path='samtools'):
    """
    Get the number of reads in a BAM file.

    Args:
        bamfile (str): Path to the BAM file.
        samtools_path (str): Path to the 'samtools' executable.

    Returns:
        int: Number of reads.
    """

    # call samtools
    out = subprocess.run([samtools_path, 'idxstats', bamfile], capture_output=True, shell=False)
    
    # read output
    out = out.stdout.decode().split('\n')
    n_reads = 0
    for l in out:
        l = l.split('\t')
        if not l == ['']:
            n_reads += int(l[2])

    return(int(n_reads/2))









