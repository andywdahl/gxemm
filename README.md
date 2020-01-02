# gxemm
R package to fit GxEMM,  a linear mixed model for polygenic GxE.

GxEMM can estimate the heritability specifc to individual environmental variables. GxEMM supports continuous and discrete environmental variables and traits. The "E" can be any context variable, including other genotypes (GxG) or inferred subtypes (for subtype-specific genetic components).

The implementation provided here will only scale to ~15,000 samples. We are currently extending the software to scale to ~50,000 samples and, in parallel, improving the underlying algorithm to scale to millions of samples.

GxEMM depends on LDAK to fit REML. We have included a version that is compatible with GxEMM in this repository. If you use this implementation of GxEMM, you must also cite the original LDAK paper:

Speed, Hemani, Johson and Balding. “Improved Heritability Estimation from Genome-wide SNPs”. The American Journal of Human Genetics 91.6 (Dec. 2012), pp. 1011–1021.

For further information on LDAK and the most current version, see http://dougspeed.com/
