# Author: cns
# KFAS backend

dlmodeler.filter.KFAS <-
		function(yt, model, raw.result=FALSE, logLik=TRUE, filter=TRUE)
{
	if(!require('KFAS')) stop("required package could not be found: KFAS")
	
	# KFAS uses ts-class representation of the time series, 
	# so each row consists of one time step, hence the transpose
	kfas.model<-KFAS::SSModel(y=t(yt),Z=model$Zt,T=model$Tt,
			R=model$Rt,H=model$Ht,Q=model$Qt,
			a1=model$a0,P1=model$P0,P1inf=model$P0inf)
	res<-KFAS::KFS(kfas.model,smoothing="none") # just filtering

	if( raw.result ) raw.res <- res else raw.res <- NA
	if( length(dim(model$Zt))==2 ) {
		# 1-step ahead prediction when observation matrix is not time-varying
		f <- kfas.model$Z[,,1] %*% res$a
	} else {
		# 1-step ahead prediction when observation matrix is time-varying
		f <- matrix(NA,NROW(res$model$Z),NCOL(yt))
		for( i in 1:NCOL(yt) ) f[,i] <- res$model$Z[,,i] %*% res$a[,i]
	}
	
	return(list(backend='KFAS',
					f=f,
					at=res$a,
					Pt=res$P,
					logLik=res$logLik,
					d=res$d,
					raw.result=raw.res))
}



dlmodeler.smooth.KFAS <-
		function(filt, raw.result=FALSE)
{
	if(!require('KFAS')) stop("required package could not be found: KFAS")
	
	res <- KFAS::KFS(filt$raw.result,smoothing="state") # same function, now smoothing
	
	if( raw.result ) raw.res <- res else raw.res <- NA
	
	return(list(backend='KFAS',
					at=res$alphahat,
					Pt=res$V,
					raw.result=raw.res))
}


