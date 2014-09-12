import sys
import pysam
from collections import Counter


# http://pysam.readthedocs.org/en/latest/api.html#pysam.AlignedRead.cigar
MATCH  = 0  # M
INS    = 1  # I
DEL    = 2  # D
SKIP   = 3  # N
SOFT   = 4  # S
HARD   = 5  # H
PAD    = 6  # P
EQUAL  = 7  # =
DIFF   = 9  # X


def cigar_profile(cigar_tuples):
	"""
	Return a dictionary that tabulates the total number
	of bases associated with each CIGAR operation.

	cigar_tuples is a list of (op, length) tuples.
	"""
	cigar_prof = Counter()
	for cigar_tuple in cigar_tuples:
		cigar_prof[cigar_tuple[0]] += cigar_tuple[1]
	return cigar_prof

def get_total_differences(cigar_prof):
	"""
	return the total number of get_total_differences
	in the alignment between the query and the reference.
	(mismatches + insertions + deletions)
	"""
	return cigar_prof[DIFF] + cigar_prof[INS] + cigar_prof[DEL]	

def get_total_unaligned(cigar_prof):
	"""
	return the total number unaligned bases (hard or softclips.)
	"""
	return cigar_prof[HARD]	+ cigar_prof[SOFT]


# header
print '\t'.join(['query', 'read_type', 'align_len', 'unalign_len', 'matches', 
	'mismatches', 'insertions', 'deletions', 'tot_errors', 'identity'])

# iterate through each BAM alignment 
# and report its alignment and error profile
bam = pysam.Samfile(sys.argv[1])
for read in bam:
	cigar_prof = cigar_profile(read.cigar)
	total_errors = get_total_differences(cigar_prof)
	unaligned_len = get_total_unaligned(cigar_prof)
	read_type = read.qname.split('_')[4]
   	print '\t'.join(str(s) for s in [read.qname, read_type, read.alen, 
   		unaligned_len, \
	   	cigar_prof[EQUAL], \
	   	cigar_prof[DIFF], \
	   	cigar_prof[INS], \
	   	cigar_prof[DEL], \
	   	total_errors, \
	   	1.0-(float(total_errors) / float(read.alen))])
