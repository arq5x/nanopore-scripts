
ONTREADS=/path/to/ONT/fast5/files/
GENOME=/path/to/reference/genome/in/fasta/format.fa
STUB="myreads"

# extract high quality reads is FASTA format.
poretools fasta --high-quality $ONTREADS > $STUB.fa

# create a LAST db of the reference genome
lastdb genome $GENOME

# align high quaklity reads to reference genome with LAST
lastal -q 1 -a 1 -b 1 genome $STUB.fa > $STUB.fa.maf

# convert the MAF to BAM with complete CIGAR (matches and mismatches)
python maf-convert.py sam $STUB.fa.maf | \
    samtools view -T $GENOME -bS - | \
    samtools sort -f - $STUB.fa.bam

# profile the errors in the alignment
python count-errors.py $STUB.fa.bam > $STUB.fa.bam.profile.txt 

# plot the error profile with R
Rscript plot-alignment-stats.R $STUB.fa.bam.profile.txt 
