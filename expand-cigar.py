import sys
import pysam
import argparse

"""Author: Aaron R. Quinlan, December 2014"""

# http://pysam.readthedocs.org/en/latest/api.html#pysam.AlignedRead.cigar
MATCH = 0  # M
INS = 1    # I
DEL = 2    # D
SKIP = 3   # N
SOFT = 4   # S
HARD = 5   # H
PAD = 6    # P
EQUAL = 7  # =
DIFF = 8   # X


def get_chrom(fasta_fh, chrom):
    """
    Return the chromosome sequence
    """
    return fasta_fh.fetch(chrom)


def expand_match(qry_seq, ref_seq):
    """
    Return an expanded list of CIGAR ops
    based upon nuceltotide matches (=, or 7)
    and mismatches (X, or 8)
    """
    prev_op = None
    curr_op = None
    length = 1
    for idx, q_nucl in enumerate(qry_seq):
        if q_nucl == ref_seq[idx]: curr_op = 7  # EQUAL (=)
        else: curr_op = 8  # DIFF (X)

        if curr_op == prev_op:
            length += 1
        elif prev_op is not None:
            yield (prev_op, length)
            length = 1
        prev_op = curr_op
    yield (curr_op, length)


def main(args):

    bam = pysam.AlignmentFile(args.bam)
    fa = pysam.FastaFile(args.fasta)
    outfile = pysam.AlignmentFile("-", "wb", template=bam)

    prev_chrom_id = None
    curr_chrom_id = None
    curr_chrom_seq = None
    for read in bam:

        # skip unaligned reads.
        if read.is_unmapped:
            continue

        curr_chrom_id = read.reference_id

        # load the current chromosome into memory
        if curr_chrom_id != prev_chrom_id:
            curr_chrom_seq = get_chrom(fa, bam.getrname(curr_chrom_id))

        # iterate through the existing CIGAR
        # and replace M ops with expanded X and = ops.
        ref_pos = read.reference_start
        qry_pos = 0
        new_cigar = []
        for cigar_tuple in read.cigar:
            op = cigar_tuple[0]
            op_len = cigar_tuple[1]
            if op == EQUAL or op == DIFF or op == MATCH:
                if op == MATCH:
                    qry_seq = read.query_sequence[qry_pos:qry_pos + op_len]
                    ref_seq = curr_chrom_seq[ref_pos:ref_pos + op_len]
                    if qry_seq == ref_seq:
                        new_cigar.append((7, len(qry_seq)))  # EQUAL (=)
                    else:  # expand the M CIGAR op into X and = ops.
                        for new_cigar_tuple in expand_match(qry_seq, ref_seq):
                            new_cigar.append(new_cigar_tuple)
                else:
                    new_cigar.append(cigar_tuple)
                ref_pos += op_len
                qry_pos += op_len
            elif op == DEL or op == SKIP:
                ref_pos += op_len
                new_cigar.append(cigar_tuple)
            elif op == INS or op == SOFT:
                qry_pos += op_len
                new_cigar.append(cigar_tuple)
            elif op == HARD:
                new_cigar.append(cigar_tuple)

        # replace the old CIGAR and write updated record.
        read.cigar=new_cigar
        outfile.write(read)
        prev_chrom_id=curr_chrom_id
    outfile.close()

if __name__ == "__main__":
    parser=argparse.ArgumentParser(prog='expand_cigar.py')
    parser.add_argument('--bam',
        dest='bam',
        metavar='STRING',
        help='The sorted BAM file whose CIGARs you wish expand.')
    parser.add_argument('--fasta',
        dest='fasta',
        metavar='STRING',
        help='The reference genome used to create the BAM.')
    args=parser.parse_args()

    main(args)

