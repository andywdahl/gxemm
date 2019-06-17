test_gxemm	<- function( K0, out_hom, out_iid, out_free, out_homE, out_homG ){

	homp	<- homp_lrt	<- iidp	<- iidp_lrt	<- h2eqp	<- sig2gp	<- sig2gp_lrt	<- sig2gp_homE	<- sig2gp_homE_lrt<- homp_freeE		<- homp_freeE_lrt<- NA

	if( !missing(out_hom) ){
	homp			<- Waldtest( out_hom$h2, out_hom$h2Covmat[1,1] )
	#homp_lrt	<- NA
	}

	if( !missing(out_hom) & !missing(out_iid) ){
	iidp			<- Waldtest( out_iid$h2[2], out_iid$h2Covmat[2,2] )
	iidp_lrt	<- LRtest( ll1=out_iid$ll, ll0=out_hom$ll, df=out_iid$df-out_hom$df )
	}

	if( !missing(out_free) ){
	h2eqp			<- h2_equal_test( out_free$h2, out_free$h2Covmat )
	sig2gp		<- MVWaldtest( out_free$sig2s[1+1:K0], out_free$sig2Var[1+1:K0,1+1:K0] )
	}

	#if( !missing(out_free) & !missing(out_iid) ){
	#sig2gp_lrt<- NA #LRtest( ll1=out_diag$ll, ll0=out_iid$ll, df=out_diag$df-out_iid$df )
	#}

	#allps	[sigtype,xval,it,'diag1']	<- h2_equal_test( out_diag$sig2s[2:4], out_diag$sig2Var[2:4,2:4] )

	if( !missing(out_homE) ){
	sig2gp_homE		<- MVWaldtest( out_homE$sig2s[1+1:K0], out_homE$sig2Var[1+1:K0,1+1:K0] )
	sig2gp_homE_lrt	<- LRtest( ll1=out_homE$ll, ll0=out_hom$ll, df=out_homE$df-out_hom$df )
	}

	if( !missing(out_homG) ){
	homp_freeE		<- Waldtest( out_homG$sig2s[1], out_homG$sig2Var[1,1] )
	homp_freeE_lrt<- LRtest( ll1=out_free$ll, ll0=out_homG$ll, df=out_free$df-out_homG$df )
	}

	list( 
		homp=homp, homp_lrt=homp_lrt,
		iidp=iidp, iidp_lrt=iidp_lrt,
		h2eqp=h2eqp, sig2gp=sig2gp, sig2gp_lrt=sig2gp_lrt,
		sig2gp_homE=sig2gp_homE, sig2gp_homE_lrt=sig2gp_homE_lrt,
		homp_freeE=homp_freeE, homp_freeE_lrt=homp_freeE_lrt
	)
}

h2_equal_test	<- function( h2s, h2Covmat ){#, prevs_sam, prevs_pop=prevs_sam ){
	K0			<- length( h2s )
	Tmat		<- cbind( 1, -diag(K0-1) )
	Tmath2	<- Tmat %*% as.numeric(h2s)
	chi2val	<- as.numeric( solve( Tmat %*% h2Covmat %*% t(Tmat) ) %*% Tmath2 ) %*% Tmath2 
	as.numeric( pchisq( chi2val, df=K0-1, lower.tail=F ) )
}

LRtest	<- function( ll1, ll0, df )
	pchisq( 2*(ll1-ll0) ,df=df , lower.tail=F )

Waldtest	<- function( z, Var )
	pnorm( z/sqrt(Var), lower.tail=F )

MVWaldtest	<- function( z, Var, eigtol=1e8 ){
	Var		<- 1/2*(Var+t(Var))
	eval	<- eigen(Var,symmetric=TRUE)$values
	if( max(eval)/(min(eval)+1e-99) > eigtol | min(eval)<0 )
		return(NA)
	pchisq( as.numeric(z) %*% (solve(Var) %*% as.numeric(z)), df=length(z), lower.tail=F )
}
