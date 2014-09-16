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



# Pass 1:
# iterate through each BAM alignment 
# and store the best alignment for the read

# best_align:
#	key is query name
#   val is tuple of (alignment length, cigar_prof)
best_align = {}
bam = pysam.Samfile(sys.argv[1])
for read in bam:
	cigar_prof = cigar_profile(read.cigar)

	if read.qname not in best_align:
		best_align[read.qname] = (read.alen, read.inferred_length, cigar_prof)
	elif read.qname in best_align \
	and read.alen > best_align[read.qname][0]:
		best_align[read.qname] = (read.alen, read.inferred_length, cigar_prof)


# Pass 2:
# Report the alignment and error profile for each read's best alignment
print '\t'.join(['query', 'read_type', 'read_len', 'align_len', 'unalign_len', 'matches', 
	'mismatches', 'insertions', 'deletions', 'tot_errors'])
for query in best_align:
	read_type = query.split('_')[4]
	alen = best_align[query][0]
	inferred_length = best_align[query][1]
	cigar_prof = best_align[query][2]
	total_errors = get_total_differences(cigar_prof)
	unaligned_len = get_total_unaligned(cigar_prof)
   	print '\t'.join(str(s) for s in [query, read_type, inferred_length, alen, 
   		unaligned_len, \
	   	cigar_prof[EQUAL], \
	   	cigar_prof[DIFF], \
	   	cigar_prof[INS], \
	   	cigar_prof[DEL], \
	   	total_errors])
