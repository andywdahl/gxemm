\name{test_gxemm}
\alias{test_gxemm}
\alias{LRtest}
\alias{Waldtest}
\alias{MVWaldtest}
\alias{h2_equal_test}

\title{Gene-x-Environnment Linear Mixed Model}
\description{Linear mixed model for polygenic gene-environment interaction}

\usage{
test_gxemm( K0, out_hom, out_iid, out_free, out_homE, out_homG )
LRtest( ll1, ll0, df )
Waldtest( z, Var )
MVWaldtest( z, Var, eigtol=1e8 )
h2_equal_test( h2s, h2Covmat )
}

\arguments{
  \item{K0}{Environment covariate matrix degrees of freedom}
  \item{out_hom,out_iid,out_free,out_homE,out_homG}{GxEMM output under different modeling assumptions}
  \item{ll1,ll0}{Alternate and null likelihoods}
  \item{df}{Degrees of freedom for the LRT}
  \item{z}{Gaussian statistic}
  \item{Var}{Variance of Gaussian statistic}
  \item{h2s,h2Covmat}{Vector of heritability estimates and their covariance matrix}
  \item{eigtol}{Maximum condition number of Var before returning NA}
}
\value{
test_gxemm returns tests comparing the various provided GxEMM models:
  \item{homp,homp_lrt}{p-values comparing Hom GxEMM model and null model}
  \item{iidp,iidp_lrt}{p-values comparing IID and Hom GxEMM models}
  \item{sig2gp,sig2gp_lrt}{p-values comparing .}
  \item{sig2gp_homE,sig2gp_homE_lrt}{p-values comparing .}
  \item{homp_freeE, homp_freeE_lrt}{p-values comparing .}
  \item{h2eqp}{p-values testing equality across discrete environmental groups under the Free GxEMM model}

LRtest, Waldtest, MVWaldtest, and h2_equal_test return p-values.
}

\author{Andy Dahl}
%\seealso{Rldak,Rgcta}
\references{Dahl, A., Cai, N. Flint, J., and Zaitlen, N. 2018. GxEMM: Extending linear mixed models to general gene-environment interactions. bioRxiv}
\details{}
\examples{

## simulate genetically homogeneous traits
Ktrue	<- 2
N	<- 2e3
L	<- 1e4
G	<- scale( matrix( rnorm(N*L), N, L ) )
K	<- 1/L*G \%*\% t(G)			# kinship matrix
Z	<- rbinom( N, 1, .2 )		# environment/subtype indicators
Z	<- cbind( Z, 1-Z )

#true h2 is .5
y	<- scale( 1/sqrt(L)*G \%*\% rnorm(L) + rnorm( N ) )

out_hom	<- GxEMM( y, X=Z[,-1], K, Z=Z, gtype='hom' )
out_iid	<- GxEMM( y, X=Z[,-1], K, Z=Z, gtype='iid' )
out_free<- GxEMM( y, X=Z[,-1], K, Z=Z, gtype='free' )

test_gxemm( K0=2, out_hom, out_iid, out_free )

}
