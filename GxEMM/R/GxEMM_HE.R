GxEMM_HE	<- function( y, X=NULL, K, Z=NULL, Ztype=c('disc','quant'), gtype=c('hom','iid','free')[1], etype=c('hom','free')[1], use_diag=TRUE ){

	Z		<- as.matrix(Z)
	X		<- as.matrix(X)
	K0	<- ncol(Z)

	X		<- cbind( 1, X )
	X_d	<- svd(X)$d
	if( max(X_d)/min(X_d) > 1e4 )
		stop( 'cbind(1,X) is nearly collinear; either X should be rescaled or a column should be dropped' )

	# residualize
	betas	<- coef( lm( y ~ -1 + X )  )
	PX	<- diag(nrow(X)) - X %*% solve( t(X)%*%X ) %*% t(X)
	y		<- as.numeric( scale( PX %*% y ) )
	y.array	<- lowtri( y %*% t(y), use_diag=use_diag )
	rm(y,X,X_d); gc()

	K.array	<- as.matrix( lowtri(K, PX, use_diag=use_diag) )
	if( gtype == 'iid' ){
		K.array		<- cbind( K.array, lowtri( K * (Z %*% t(Z)), PX, use_diag=use_diag ) )
	} else if( gtype == 'free' ){
		for( k in 1:K0 )
			K.array	<- cbind( K.array, lowtri( K * (Z[,k,drop=F] %*% t(Z[,k,drop=F])), PX, use_diag=use_diag ) )
	}
	rm(K); gc()

	if( etype == 'iid' & Ztype == 'disc' ){
		K.array		<- cbind( K.array, lowtri( diag(diag(Z %*% t(Z))), PX, use_diag=use_diag ) )
	} else if( etype == 'free' )
		for( k in 1:ifelse( Ztype == 'disc', K0-1, K0 ) )
			K.array	<- cbind( K.array, lowtri( diag( as.numeric( Z[,k]^2 ) ), PX, use_diag=use_diag ) )
	K.array			<- cbind( K.array, lowtri(diag(nrow(K)), PX, use_diag=use_diag) )
	rm( PX ); gc()

	sig2s		<- solve( t(K.array) %*% K.array ) %*% ( t(K.array) %*% y.array )

	sig2out	<- sig2map_HE( sig2s, gtype, etype, K0, disc_Z=disc_Z )

	list( h2=sig2out$h2, sig2g=sig2out$sig2g, sig2e=sig2out$sig2e, df=ncol(K.array), sig2s=sig2s, ll=NA, h2Covmat=NA, sig2Var=NA, betas=betas )
}

lowtri	<- function(V,PX,use_diag=TRUE){
	if( ! missing(PX) )
		V	<- PX %*% V %*% PX
	if( use_diag ){
		return(c(V[lower.tri(V,diag=FALSE)],diag(V)/2))
	} else {
		return(  V[lower.tri(V,diag=FALSE)])
	}
}
