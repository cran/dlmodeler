# Author: cns
# dlm backend

dlmodeler.filter.dlm <-
		function(yt, model, raw.result=FALSE, logLik=FALSE, filter=TRUE)
{
	if(!require('dlm')) stop("required package could not be found: dlm")
	Ttvar <- length(dim(model$Tt))==3
	Rtvar <- length(dim(model$Rt))==3
	Qtvar <- length(dim(model$Qt))==3
	Ztvar <- length(dim(model$Zt))==3
	Htvar <- length(dim(model$Ht))==3
	if( Ttvar|Rtvar|Qtvar|Htvar ) stop("This kind of time varying model is currently unsupported")
	# TODO handle more time varying cases
	Ht <- model$Rt %*% model$Qt %*% t(model$Rt)
	if( Ztvar ) {
		n1 <- NROW(model$Zt)
		n2 <- NCOL(model$Zt)
		n3 <- dim(model$Zt)[3]
		FF <- matrix(0,n1,n2)
		JFF <- matrix(1:(n1*n2),n1,n2)
		X <- matrix(NA,n3,n1*n2)
		for( i in 1:n1 ) for( j in 1:n2 ) X[,i+(j-1)*n1] <- model$Zt[i,j,]
		mdlm <- dlm::dlm(
				m0=model$a0,
				C0=model$P0+1e7*model$P0inf,
				FF=FF,JFF=JFF,X=X,
				V=model$Ht,
				GG=model$Tt,
				W=Ht)
	} else {
		mdlm <- dlm::dlm(
				m0=model$a0,
				C0=model$P0+1e7*model$P0inf,
				FF=model$Zt,
				V=model$Ht,
				GG=model$Tt,
				W=Ht)
	}
	# dlm has currently no way of computing logLik and filtering at the
	# same time, so dlmFilter and dlmLL are only called on request
	if( filter ) {
		res <- dlm::dlmFilter(t(yt),mdlm,simplify=!raw.result)
		res.Pt <- dlm::dlmSvd2var(res$U.R,res$D.R)
		Pt <- array(NA,dim=c(NROW(Ht),NCOL(Ht),length(res.Pt)))
		for( i in 1:length(res.Pt) ) Pt[,,i] <- res.Pt[[i]]
	} else {
		res <- NA
		Pt <- NA
	}
	if( raw.result ) raw.res <- res else raw.res <- NA
	if( logLik ) logLik <- dlm::dlmLL(yt,mdlm) else logLik <- NA
	
	## if( length(dim(model$Zt))==2 ) {
	##     # 1-step ahead prediction when observation matrix is not time-varying
	##     fm <- model$Zt %*% t(res$m)
	##     fa <- model$Zt %*% t(res$a)
	## } else {
	##     # 1-step ahead prediction when observation matrix is time-varying
	##     fm <- matrix(NA,NROW(model$Zt),NCOL(yt))
	##     fa <- matrix(NA,NROW(model$Zt),NCOL(yt))
	##     for( i in 1:NCOL(yt) ) fm[,i] <- model$Zt[,,i] %*% t(res$m[,i])
	##     for( i in 1:NCOL(yt) ) fa[,i] <- model$Zt[,,i] %*% t(res$a[,i])
	## }
	
	return(list(backend='dlm',
					f=t(res$f), ## fm=fm, fa=fa,
					at=t(res$a),
					Pt=Pt,
					logLik=logLik,
					d=0,
					raw.result=raw.res))
}



dlmodeler.smooth.dlm <-
		function(filt, raw.result=FALSE)
{
	if(!require('dlm')) stop("required package could not be found: dlm")
	
	res <- dlm::dlmSmooth.dlmFiltered(filt$raw.result)
	
	if( raw.result ) raw.res <- res else raw.res <- NA
	res.Pt <- dlm::dlmSvd2var(res$U.S,res$D.S)
	Pt <- array(NA,dim=c(NCOL(res$D.S),NCOL(res$D.S),length(res.Pt)))
	for( i in 1:length(res.Pt) ) Pt[,,i] <- res.Pt[[i]]
	
	return(list(backend='dlm',
					at=t(res$s),
					Pt=Pt,
					raw.result=raw.res))
}


