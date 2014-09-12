library('ggplot2')

mapstat <- read.table('reads.fa.bam.profile.txt', header=TRUE)

ggplot(mapstat, 
       aes(x=identity, 
           y=align_len/(align_len+unalign_len), 
           colour=read_type)) +
  xlab("Percent identity (mm+ins+del)/(align_len)") +
  ylab("Fraction of read aligned") +
  geom_point(alpha=0.8) + 
  xlim(0.4,1.0) +
  theme_bw(base_size=22) 
