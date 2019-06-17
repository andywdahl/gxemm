GxEMM	<- function(
	y, X=NULL, K, Z=NULL, # if Z is NULL, just fits GREML
	gtype=c('hom','iid','free')[1], etype=c('hom','free')[1],
	binary=FALSE, prev,
	ldak_loc='~/GxEMM/code/ldak5.linux ', lowmem=FALSE,
	tmpdir=paste0( 'gxemm_tmp_', round(runif(1)*1e5) ), # temporary kinship matrix files to feed to LDAK for all the K * ZZT terms (or I * ZZT for the heterogeneous noise option, which GCTA does not do)
	keep_ldak=FALSE # helpful for debugging failed GxEMM runs
){
	Z		<- as.matrix(Z)
	X		<- as.matrix(X)
	K0	<- ncol(Z)

	#### check not already ran
	if( file.exists( paste0( tmpdir, '/y.phen' ) ) ) stop('Output file already exists')

	#### set up LDAK output
	system( paste0( 'mkdir ', tmpdir ) )
	if(keep_ldak){
		outdir=paste0( 'gxemm_out_', round( abs(rnorm(1) + X[sample(10,1)] + y[sample(10,1)] )*1e5, 0 ) )
		system( paste0( 'mkdir ', outdir ) )
	} else {
		outdir=tmpdir
	}

	### write pheno + covars
	if( binary ){
		stopifnot( all( is.na(y) | (y%in%1:2) ) )
		if( missing( prev ) ){
			warning( 'Prevalence missing. Assuming binary trait isn\'t scertained, which fails badly for rare disease case/control data' )
			prev	<- sum(y==2)/sum(y%in%1:2)
		}
	} else {
		y	<- scale(y)
	}
	write_pheno( y, file=paste0( tmpdir, '/y.phen')	 )
	write_pheno( X, file=paste0( tmpdir, '/X.qcovar') )
	X			<- cbind( 1, X )

	X_singvals	<- svd(X)$d
	if( max(X_singvals)/min(X_singvals) > 1e4 )
		stop( 'cbind(1,X) is nearly collinear; either X should be scaled or a column should be dropped' )

	### write each of the kinship matrices
	ws	<- write_kin( tmpdir, K, index=1, ldak_loc, X ) # hom K
	if( gtype == 'iid' ){
		ws	<- c( ws, write_kin( tmpdir, K * (Z %*% t(Z)), index=2, ldak_loc, X ) ) # iid K
	} else if( gtype == 'free' ){
		for( k in 1:K0 )
			ws	<- c( ws, write_kin( tmpdir, K * ( Z[,k,drop=F] %*% t(Z[,k,drop=F])	), index=1+k, ldak_loc, X ) ) # K * (Zk Zk^T)
	}
	rm(K); gc()
	if( etype == 'free' )
		for( k in 1:ifelse( binary, K0, K0-1 ) ) #K0-th is omitted because singular with Hom noise for categorical Z. In simulations with binary traits, however, the background seems to play a different role and remains identified for categorical Z, so I use all K0 K*E matrices in the Free mode. This is very heuristic at this stage, but the simulation results are consistent across many simulation settings.
			ws	<- c( ws, write_kin( tmpdir, diag( Z[,k]^2 ), index=length(ws)+1, ldak_loc, X ) ) # I * (Zk Zk^T)
	rm(Z); gc()
	r		<- length(ws)
	ws	<- c( ws, 1 - ncol(X)/nrow(X) )
	write.table( paste0( tmpdir, '/K.', 1:r ), file=paste0( tmpdir, '/multi_grm.txt' ), sep='\t', row.names=F, col.names=F, quote=F )

	### adjust kinships for the covariates. I'm not sure this is exactly the right way to adjust for covariates, but it's at least close
	if( binary )
		for( i in 1:r )
			adjust_ldak_kinship( prefix=paste0( tmpdir, '/K.', i ), ldak_loc, covarfile=paste0( tmpdir, '/X.qcovar ' ) )

	failed	<- TRUE
	try({
		# R -> LDAK
		system( paste0( ldak_loc,
			ifelse( r==1, 
				paste0(	' --grm ', tmpdir, '/K.1 '  ),
				paste0( ' --mgrm ', tmpdir, '/multi_grm.txt ' )
			), 
			' --pheno ', tmpdir, '/y.phen ',
			' --covar ', tmpdir, '/X.qcovar ',
			ifelse( binary, ' --pcgc ', ' --reml ' ), outdir, '/tmp',
			ifelse( binary, paste0( ' --prevalence ', prev ), '' ),
			ifelse( binary, '', ' --reml-iter 1000 --tolerance 1e-6' ),
			' --memory-save ', ifelse( lowmem , 'YES ', 'NO ' ),
			' --kinship-details NO'
		))
		# LDAK -> R
		gout	<- extract_ldak( binary, outdir, ws, r )
		failed	<- FALSE
	})
	system( paste0( 'rm -rf ', tmpdir ) )
	if( failed ) return('FAILED')
	if( ! binary & gout$niter == 1000 ){ print('LDAK did not converge'); return('LDAK did not converge') }

	# convert to more meaningful parameters
	sig2out		<- sig2map( gout$sig2s, gtype, etype, K0, binary=binary )
	h2Covmat	<- msm::deltamethod( sig2out$form, gout$sig2s, gout$sig2Var, ses=F )

	list( h2=sig2out$h2, sig2g=sig2out$sig2g, sig2e=sig2out$sig2e, df=r, sig2s=gout$sig2s, ll=gout$ll, h2Covmat=h2Covmat, sig2Var=gout$sig2Var, betas=gout$betas )
}
