library('ggplot2')

args <- commandArgs(trailingOnly = TRUE)

mapstat <- read.table(args[1], header=TRUE)

ggplot(mapstat, 
       aes(x=100*matches/(matches+deletions+insertions+mismatches), 
           y=(align_len+insertions-deletions)/(align_len+unalign_len-deletions+insertions), 
           colour=read_type)) +
  xlab("Percent identity: 100* matches/(matches+deletions+insertions+mismatches)") +
  ylab("Fraction of read aligned: \n(align_len+ins-del)/(align_len+unalign_len-del+ins)") +
  geom_point(alpha=0.5) + 
  xlim(20,100) +
  ylim(0.0,1.0) +
  theme_bw(base_size=24)