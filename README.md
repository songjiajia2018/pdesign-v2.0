# pdesign-v2.0
a new version of design
Since probes with high stability of secondary structures are difficult to hybridize with the 16S sequences, candidate probes should avoid with low free energy. In this study, we further optimized the pdesign, an algorithm developed in our previous work, to prevent the formation of secondary DNA structures.
a k-mer-based algorithm
first use:
gzip -d spool.txt.gz

usage:
perl pdesign.pl [order name] [target name] [output former] [thread] >[final output name]

example:
perl pdesign.pl Porphyromonadaceae Barnesiella barn 10 >barn-result
