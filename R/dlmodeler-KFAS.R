# Author: cns
# KFAS backend

dlmodeler.filter.KFAS <-
		function(yt, model, raw.result=FALSE, logLik=TRUE, filter=TRUE)
{
	if(!require('KFAS')) stop("required package could not be found: KFAS")
	res <- KFAS::kf(yt=yt,Zt=model$Zt,Tt=model$Tt,
			Rt=model$Rt,Ht=model$Ht,Qt=model$Qt,
			a1=model$a0,P1=model$P0,P1inf=model$P0inf,
			optcal=c(FALSE,FALSE,FALSE,FALSE))
	if( raw.result ) raw.res <- res else raw.res <- NA
	if( length(dim(model$Zt))==2 ) {
		# 1-step ahead prediction when observation matrix is not time-varying
		f <- res$Zt %*% res$at
	} else {
		# 1-step ahead prediction when observation matrix is time-varying
		f <- matrix(NA,NROW(res$Zt),NCOL(yt))
		for( i in 1:NCOL(yt) ) f[,i] <- res$Zt[,,i] %*% res$at[,i]
	}
	return(list(backend='KFAS',
					f=f,
					at=res$at,
					Pt=res$Pt,
					logLik=res$lik,
					d=res$d,
					raw.result=raw.res))
}



dlmodeler.smooth.KFAS <-
		function(filt, raw.result=FALSE)
{
	if(!require('KFAS')) stop("required package could not be found: KFAS")
	res <- KFAS::ks(filt$raw.result)
	if( raw.result ) raw.res <- res else raw.res <- NA
	return(list(backend='KFAS',
					at=res$ahat,
					Pt=res$Vt,
					raw.result=raw.res))
}


