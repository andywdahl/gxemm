sig2map	<- function( sig2s, gtype, etype, K0, binary=FALSE ){
	if( gtype == 'hom' & etype == 'hom' ){
		sig2g	<- sig2s[1]
		sig2e	<- sig2s[2]
		form	<- list( as.formula("~(x1)/(x1+x2)"))

	} else if( gtype == 'iid' ){
		if( etype == 'free' ) stop('iid + het noise is not very meaningful')
		sig2g	<- sig2s[1:2]
		sig2e	<- sig2s[3]
		form	<- lapply( 1:2, function(k) as.formula( paste0( '~x', k, '/(x1+x2+x3)' ) ) )

	} else if( gtype == 'hom' & etype == 'free' ){
		sig2g		<- sig2s[1]
		if( binary ){
			sig2e	<- sig2s[1+1:K0]+sig2s[K0+2]
			form	<- lapply( 1:K0, function(k) as.formula( paste0( '~x1/(x1+x',1+k,'+x',K0+2,')' ) ) )
		} else {
			sig2e	<- c(sig2s[1+1:(K0-1)],0)+sig2s[1+K0]
			form	<- lapply( 1:K0, function(k) as.formula( paste0( '~x1/(x1+x',ifelse( k==K0, paste0(1+k,'+x'), '' ),K0+1,')') ) )
		}

	} else if( gtype == 'free' & etype == 'free' ){
		sig2g		<- sig2s[1:(K0+1)]
		if( binary ){
			sig2e	<- sig2s[K0+1+1:K0]+sig2s[2*K0+2]
			form	<- lapply( 1:K0, function(k) as.formula( paste0( '~(x1+x',1+k,')/(x1+x',1+k,'+x',1+K0+k,'+x',2*K0+2,')') ) )
		} else {
			sig2e	<- c(sig2s[K0+1+1:(K0-1)],0)+sig2s[2*K0+1]
			form	<- lapply( 1:K0, function(k) as.formula( paste0( '~(x1+x',1+k,')/(x1+x',1+k,'+x',ifelse(k==K0,'',paste0(1+K0+k,'+x')),2*K0+1,')') ) )
		}

	} else if( gtype == 'free' & etype == 'hom' ){
		sig2g	<- sig2s[1:(K0+1)]
		sig2e	<- sig2s[1+K0+1]
		form	<- lapply( 1:K0, function(k) as.formula(paste0( '~(x1+x', 1+k, ')/(x1+x', 1+k, '+x', 1+K0+1, ')')) )

	} else {
		stop(c(gtype,etype))
	}

	if( gtype == 'hom' & etype == 'hom' ){
		h2				<- sig2g/(sig2g+sig2e)
		names(h2)	<- 'hom'
	} else if( gtype == 'iid' ){
		h2				<- sig2g/(sum(sig2g)+sig2e)
		names(h2)	<- c( 'hom', 'het' )
	} else {
		if( gtype == 'free' ){
			h2					<- (sig2g[1] + sig2g[-1])/(sig2g[1] + sig2g[-1] + sig2e)
			names(sig2g)<- c( 'sig2g_hom', paste0( 'sig2g_', 1:K0 ) )
		} else {
			h2					<- sig2g/(sig2g+sig2e)
			names(sig2g)<- 'sig2g_hom'
		}
		names(h2)	<- paste0( 'h2_', 1:K0 )
	}
	list( sig2g=sig2g, sig2e=sig2e, form=form, h2=h2	)
}
