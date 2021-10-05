extract_ldak	<- function( pcgc, outdir, ws, r ){
	if( ! pcgc ){

		out			<- read.table( paste0( outdir, '/tmp.vars' ), header=TRUE ) 
		sig2s		<- as.numeric( as.character( out[,'Variance'] ) )
		sig2ses	<- as.numeric( as.character( out[,'SD'] ) )
		rm( out )

		out			<- as.matrix( read.table( paste0( outdir, '/tmp.cross' ), fill=T, header=TRUE ) )
		sig2Var	<- rbind( cbind( out, -rowSums(out) ), c( -rowSums(out), sum(out) ) )
		sig2Var	<- 1/2*( sig2Var + t(sig2Var) )
		rm( out )

		### scale back up
		out		<- read.table( paste0( outdir, '/tmp.reml' ), fill=T ) 
		ll		<- as.numeric( as.character( out[which(out[,1]=='Alt_Likelihood')	,2] ) )
		ll0		<- as.numeric( as.character( out[which(out[,1]=='Null_Likelihood'),2] ) )
		rm( out )

		out		<- read.table( paste0( outdir, '/tmp.coeff' ), fill=T, header=T ) 
		betas			<- as.numeric(out[,'Effect'])
		beta.ses	<- as.numeric(out[,'SD'])
		beta.ps		<- as.numeric(out[,'P'])
		rm( out )

		niter	<- nrow( read.table( paste0( outdir, '/tmp.progress' ) )) - 2

	} else {

		niter	<- NA
		ll		<- NA
		#TODO:
		betas			<- NA
		beta.ses	<- NA

		out			<- read.table( paste0( outdir, '/tmp.pcgc' ), fill=T, col.names=1:6 )
		sig2s		<- as.numeric( as.character( out[which( out[,1] %in% paste0( 'Her_K', 1:r ) ),2] ) )
		sig2ses	<- as.numeric( as.character( out[which( out[,1] %in% paste0( 'Her_K', 1:r ) ),3] ) )
		rm( out )

		sig2s		<- c( sig2s, 1-sum(sig2s) )
		sig2ses	<- c( sig2ses, 0 )

		#covh2	<- as.numeric( as.character( out[which( out[,1] == 'Covar_Heritability' ),2] ) )
		sig2Var		<- diag( sig2ses^2 ) ### TODO: read in full Covariance matrix from new version of LDAK!

	}

	### rescaling
	sig2s			<- sig2s	 / ws
	sig2ses		<- sig2ses / ws
	sig2Var	<- diag(1/ws) %*% sig2Var %*% diag(1/ws)

	return( list( sig2s=sig2s, sig2ses=sig2ses, ll=ll, sig2Var=sig2Var, beta.ses=beta.ses, betas=betas, niter=niter ) )
}

adjust_ldak_kinship	<- function(prefix,ldak_loc,covarfile)
	system( paste0( ldak_loc, ' --adjust-grm ', prefix, ' --grm ', prefix, ' --covar ', covarfile, ' --kinship-details NO' ) )

write_pheno	<- function(y,file,IDs=1:nrow(as.matrix(y)))
	write.table( cbind(IDs,IDs,y), file=file, sep='\t', row.names=F, col.names=F, quote=F )

write_kin	<- function(tmpdir,K,index,ldak_loc, X0, IDs=1:nrow(K) ){

	# probably unnecessary if covariates are handled elsewhere?
	Xa		<- X0 %*% solve(t(X0)%*%X0)
	KXXa	<- (K %*% X0) %*% t(Xa)
	K			<- K - KXXa - t(KXXa) + X0 %*% ( t(Xa) %*% ( KXXa ) )
	rm( X0, Xa, KXXa )

	N	<- nrow(K)
	w	<- mean(diag(K))
	K	<- K/w
	prefix	<- paste0( tmpdir, '/K.', index )
	keep = which( lower.tri(K,diag=T) , arr.ind=T )
	non_mis = rep( 1  , length(keep[,1]) )
	keep = cbind( keep , non_mis , K[ keep ] )
	rm( K ); gc()
	keep = keep[ order(keep[,1],keep[,2]) , ]
	write.table( keep						, file=gzfile(paste0( prefix, '_tmp.grm.gz' )), col.names=F, row.names=F, quote=F )
	write.table( cbind(1:N,1:N)	, file=       paste0( prefix, '_tmp.grm.id' ) , col.names=F, row.names=F, quote=F )
	system( paste0( ldak_loc, ' --convert-gz ', prefix, ' --grm ', paste0(prefix,'_tmp') ) )
	system( paste0( 'rm -rf ',																						prefix,'_tmp.*') )
	w
}
