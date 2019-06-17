\name{GxEMM}
\alias{GxEMM}
\alias{GxEMM_HE}

\title{Gene-x-Environnment Linear Mixed Model}
\description{Linear mixed model for polygenic gene-environment interaction}

\usage{
GxEMM(
	y, X=NULL, K, Z=NULL,
	gtype=c('hom','iid','free')[1], etype=c('hom','free')[1],
	binary=FALSE, prev,
	ldak_loc='~/GxEMM/code/ldak5.linux ', lowmem=FALSE,
	tmpdir=paste0( 'gxemm_tmp_', round(runif(1)*1e5) ),
	keep_ldak=FALSE
)
GxEMM_HE(
	y, X=NULL, K, Z=NULL,
	gtype=c('hom','iid','free')[1], etype=c('hom','free')[1], use_diag=TRUE
)
}

\arguments{
  \item{y}{A quantitative trait vector}
  \item{X}{A matrix of homogeneous covariates}
  \item{K}{A Kinship matrix}
  \item{Z}{A matrix of environmental covariates}
  \item{gtype}{Controls types of gene-environment interactions assumed}
  \item{etype}{Controls types of noise-environment interactions assumed}
  \item{binary}{Whether y is binary or quantitative}
  \item{prev}{For binary traits, the population trait prevalence}
  \item{ldak_loc}{Path to LDAK code}
  \item{lowmem}{Whether to run --save-mem YES in LDAK}
  \item{tmpdir}{Temporary directory to write e.g. GRMs}
  \item{keep_ldak}{Whether to remove the folder containing LDAK output}
  \item{use_diag}{Whether to include diagonal entries of phenotypic covariance matrix in HE objective function}
}
\value{
  \item{h2}{Heritability estimates}
  \item{sig2g}{Genetic variance estimates}
  \item{sig2e}{Noise variance estimates}
  \item{df}{Number of fitted variance components}
  \item{sig2s}{Variance component estimates}
  \item{ll}{Log-likelihood for the fitted gxemm}
  \item{h2Covmat}{Covariance matrix for heritability estimates}
  \item{sig2Var}{Variance component estimator covariance matrix}
  \item{betas}{Fixed effect estimates for X}
}

\author{Andy Dahl}
\references{Dahl, A., Cai, N. Flint, J., and Zaitlen, N. 2018. GxEMM: Extending linear mixed models to general gene-environment interactions. bioRxiv}
\note{X must not include an intercept, in the sense that cbind( 1, X ) must have full rank. If GxE components are fitted, X should include columns of Z to adjust for main environmental effects. y is internally scaled to mean 0, variance 1.}
\details{Fits a GxEMM to trait y with covariates X, relatedness matrix K and environments Z. `gtype'/`etype' determine which genetic/noise variance components are fit: either only the homogeneous effect (the ordinary kinship for `gtype', the identity matrix for `etype'); or adding a single GxE kinship matrix assuming equal heritability explained by each environment (`gtype=iid'); or adding a free environment-specfic contributions (`free'). See Dahl et al 2018 for further details.

The output takes different forms under different `gtype'/`etype' choices.

If both are `hom', then h2, sig2g, and sig2e are each scalar parameters.

If `gtype' is `iid' and `etype' is `hom', the h2 and sig2g parameters have length 2, for the homogeneous/heterogeneous genetic variance components, respectively.

If `gtype' is `free', then sig2g has length `1+ncol(Z)', with first entry giving the homogeneous genetic variance and the remainder giving the heterogeneous heritability specific to environment interactions with each column of `Z'--these latter terms are only easily interpretable for discrete environments.

If `etype' is `free', then sig2e has length `ncol(Z)', with each entry giving the noise level in each environment indicated by `Z'. (For non-discrete environments, the interpretation is subtler and the homogeneous noise level is identified; see Dahl et al 2018 for details.)

Finally, if either `gtype' or `etype' is `free', then h2 is length `1+ncol(Z)', analogous to sig2g, with first entry giving the homogeneous heritability and the rest giving the heterogeneous heritability specific to environment interactions with each column of `Z'.
}

\examples{

## simulate genetically homogeneous traits
Ktrue	<- 2
N	<- 2e3
L	<- 1e3
G	<- scale( matrix( rnorm(N*L), N, L ) )
K	<- 1/L*G \%*\% t(G)			# kinship matrix
Z	<- rbinom( N, 1, .2 )		# environment/subtype indicators
Z	<- cbind( Z, 1-Z )

#true h2 is .5
y	<- scale( 1/sqrt(L)*G \%*\% rnorm(L) + rnorm( N ) )

####################################
#### homogeneous/ordinary GCTA ####
####################################
out_hom	<- GxEMM( y, X=Z[,-1], K, Z=Z, gtype='hom' )
out_hom$h2

out_hom_HE	<- GxEMM_HE( y, X=Z[,-1], K, Z=Z, gtype='hom' )
out_hom_HE$h2

####################################
#### simple GxE model (slight generalization of GxE model in Yang 2010/GCTA)
####################################
out_iid	<- GxEMM( y, X=Z[,-1], K, Z=Z, gtype='iid' )
out_iid$h2

### compare hom to simple model:
lrt     <-	out_iid$ll - out_hom$ll
.5*pchisq( 2*lrt, df=1, lower.tail=FALSE )	 # Self + Liang 1980-ish

####################################
#### free GxE model ####
####################################
out_free<- GxEMM( y, X=Z[,-1], K, Z=Z, gtype='free' )
out_free$h2

### compare hom to full model:
lr	<- out_free$ll - out_hom$ll
pchisq( 2*lr, df=(Ktrue-1), lower.tail=FALSE )

### compare simple to full model:
lr	<- out_free$ll - out_iid$ll
pchisq( 2*lr, df=(Ktrue-1), lower.tail=FALSE )

out_free<- GxEMM( y, X=Z[,-1], K, Z=Z, gtype='free', etype='free' )
out_free$h2
}
