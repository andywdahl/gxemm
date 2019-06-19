# gxemm
Linear mixed model to fit polygenic GxE.

GxEMM can estimate the heritability specifc to individual environmental variables. GxEMM supports continuous and discrete environmental variables and traits.

The implementation provided here will only scale to ~15,000 samples.

GxEMM depends on LDAK to fit REML. We have included a version that is compatible with GxEMM in this repository. If you use this implementation of GxEMM, you must also cite the original LDAK paper:

Speed, Hemani, Johson and Balding. “Improved Heritability Estimation from Genome-wide SNPs”. The American Journal of Human Genetics 91.6 (Dec. 2012), pp. 1011–1021.

For further information on LDAK and the most up-to-date version, see: http://dougspeed.com/
