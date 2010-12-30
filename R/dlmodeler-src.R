
# Internal list of things to do:
#- TODO: sspir support
#- TODO: spline
#- TODO: SARIMA
#- TODO: generalized dlm support
#- TODO: plot functions and diagnostic checking
#- TODO: automatic forecasting
#- TODO: formulas (like package sspir)
#- TODO: exponential smoothing
#- TODO: vignette

 ###################
# general functions #
 ###################

dlmodeler.build <- 
function(a0=NULL, P0=NULL, P0inf=NULL,
         Tt=NULL, Rt=NULL, Qt=NULL,
         Zt=NULL, Ht=NULL,
         dimensions=NULL,
         name='noname', components=list())
{
	if( is.null(dimensions) )
	{
		components[[name]] <- matrix(1,NROW(Tt),1)
		mod <- list(
			a0=a0, P0=P0, P0inf=P0inf,
			Tt=Tt, Rt=Rt, Qt=Qt,
			Zt=Zt, Ht=Ht,
			name=name,
			components=components)
	}
	else
	{
		if( length(dimensions)!=3 ) stop("dimensions are invalid, they should have the form: c(m,r,d)")
		m <- dimensions[1] # dim state vector
		r <- dimensions[2] # dim state disturbances
		d <- dimensions[3] # dim observation vector
		components[[name]] <- matrix(1,m,1)
		mod <- list(
			a0=matrix(0,m,1),
			P0=matrix(0,m,m),
			P0inf=matrix(0,m,m),
			Tt=matrix(0,m,m),
			Rt=matrix(0,m,r),
			Qt=matrix(0,r,r),
			Zt=matrix(0,d,m),
			Ht=matrix(0,d,d),
			name=name,
			components=components)
	}
	class(mod) <- 'dlmodeler' # a nice class name!
	return(mod)
}



dlmodeler.check <-
function(model, yt=NULL)
{
	if( class(model)!='dlmodeler' ) stop("model should be of class 'dlmodeler'")
	t <- FALSE        # time varying model
	n <- rep(NA,5)    # number of time steps (in case of time-varying model)
	m <- NROW(model$Tt) # dim of state vector
	r <- NROW(model$Qt) # dim of state disturbance
	d <- NROW(model$Zt) # dim of observation vector
	s <- TRUE         # status (TRUE if model is valid)
	# check the dimensions of the initial conditions
	if( NROW(model$a0)!=m ) { warning("Inconsistent model: dimension of a0 should be ",m," x 1"); s <- FALSE }
	if( NCOL(model$a0)!=1 ) { warning("Inconsistent model: dimension of a0 should be ",m," x 1"); s <- FALSE }
	if( NROW(model$P0)!=m ) { warning("Inconsistent model: dimension of P0 should be: ",m," x ",m); s <- FALSE }
	if( NCOL(model$P0)!=m ) { warning("Inconsistent model: dimension of P0 should be: ",m," x ",m); s <- FALSE }
	if( NROW(model$P0inf)!=m ) { warning("Inconsistent model: dimension of P0inf should be: ",m," x ",m); s <- FALSE }
	if( NCOL(model$P0inf)!=m ) { warning("Inconsistent model: dimension of P0inf should be: ",m," x ",m); s <- FALSE }
	# check the dimensions of the state equation
	if( NROW(model$Tt)!=m ) { warning("Inconsistent model: dimension of Tt should be: ",m," x ",m); s <- FALSE }
	if( NCOL(model$Tt)!=m ) { warning("Inconsistent model: dimension of Tt should be: ",m," x ",m); s <- FALSE }
	if( NROW(model$Rt)!=m ) { warning("Inconsistent model: dimension of Rt should be: ",m," x ",r); s <- FALSE }
	if( NCOL(model$Rt)!=r ) { warning("Inconsistent model: dimension of Rt should be: ",m," x ",r); s <- FALSE }
	if( NROW(model$Qt)!=r ) { warning("Inconsistent model: dimension of Qt should be: ",r," x ",r); s <- FALSE }
	if( NCOL(model$Qt)!=r ) { warning("Inconsistent model: dimension of Qt should be: ",r," x ",r); s <- FALSE }	
	# check the dimensions of the observation equation
	if( NROW(model$Zt)!=d ) { warning("Inconsistent model: dimension of Zt should be: ",d," x ",m); s <- FALSE }
	if( NCOL(model$Zt)!=m ) { warning("Inconsistent model: dimension of Zt should be: ",d," x ",m); s <- FALSE }
	if( NROW(model$Ht)!=d ) { warning("Inconsistent model: dimension of Ht should be: ",d," x ",d); s <- FALSE }
	if( NCOL(model$Ht)!=d ) { warning("Inconsistent model: dimension of Ht should be: ",d," x ",d); s <- FALSE }
	# check the time-varying components
	if( length(dim(model$Tt))==3 ) { n[1] <- dim(model$Tt)[3]; t <- TRUE }
	if( length(dim(model$Rt))==3 ) { n[2] <- dim(model$Rt)[3]; t <- TRUE }
	if( length(dim(model$Qt))==3 ) { n[3] <- dim(model$Qt)[3]; t <- TRUE }
	if( length(dim(model$Zt))==3 ) { n[4] <- dim(model$Zt)[3]; t <- TRUE }
	if( length(dim(model$Ht))==3 ) { n[5] <- dim(model$Ht)[3]; t <- TRUE }
	if( t ) if (min(n,na.rm=TRUE)!=max(n,na.rm=TRUE)) { warning("Inconsistent model: dimension of time-varying terms (Tt, Rt, Qt, Zt, Ht) should match: ",n); s <- FALSE }
	# check the components
	for( compname in names(model$components) ) {
		comp <- model$components[[compname]]
		if( NROW(comp)!=m ) { warning("Inconsistent model: dimension of component ",compname," should be ",m," x 1"); s <- FALSE }
		if( NCOL(comp)!=1 ) { warning("Inconsistent model: dimension of component ",compname," should be ",m," x 1"); s <- FALSE }
		if( apply(comp==1|comp==0,2,"sum")!=m ) { warning("Inconsistent model: ",compname," be a vector of zeros and ones"); s <- FALSE }
	}
	if( apply(model$components[[model$name]]==1,2,"sum")!=m ) { warning("Inconsistent model: component named ",model$name," should be a vector of ones"); s <- FALSE }
	# check compatibility with data
	if( !is.null(yt) ) {
		if( !is.matrix(yt) ) yt <- matrix(yt,nrow=1)
		if( NROW(yt)!=d ) { warning("Inconsistent data: dimension of yt should be: ",d," x ",'n'); s <- FALSE }
		if( t & NCOL(yt)!=max(n,na.rm=TRUE) ) { warning("Inconsistent data: dimension of yt should be: ",d,"x",max(n,na.rm=TRUE)); s <- FALSE }
	}
	return(c(s,m,r,d,t,n))
}



print.dlmodeler <-
function(x,...)
{
	mdl.dim <- dlmodeler.check(x)
	if( mdl.dim[5]==FALSE ) {
		cat("constant dlmodel(state dim=",mdl.dim[2],", dist dim=",mdl.dim[3],", obs dim=",mdl.dim[4],") '",x$name,"'\n",sep='')
	} else {
		n <- max(mdl.dim[6:10],na.rm=TRUE)
		cat("time-varying dlmodel(state dim=",mdl.dim[2],", dist dim=",mdl.dim[3],", obs dim=",mdl.dim[4],", time steps=",n,") '",x$name,"'\n",sep='')
		if( !is.na(mdl.dim[6]) ) cat(" - state transition matrix is time-variant\n")
		if( !is.na(mdl.dim[7]) ) cat(" - state disturbance selection matrix is time-variant\n")
		if( !is.na(mdl.dim[8]) ) cat(" - state disturbance covariance matrix is time-variant\n")
		if( !is.na(mdl.dim[9]) ) cat(" - observation design matrix is time-variant\n")
		if( !is.na(mdl.dim[10]) ) cat(" - observation disturbance covariance matrix is time-variant\n")
	}
	nb.components <- length(x$components)
	if( nb.components>1 ) {
		cat(" - model has",nb.components,"components: ")
		cat(names(x$components),sep=", ")
		cat("\n")
	} else if( nb.components==1 ) {
		cat(" - model has",nb.components,"component: ")
		cat(names(x$components),sep=", ")
		cat("\n")
	}
	invisible(x)
}



dlmodeler.extract <-
function(fs, model, compnames=NULL, type=c("observation","state"), value=c("mean","covariance","interval"), prob=.90)
{
	if( class(fs)!='dlmodeler.filtered' & class(fs)!='dlmodeler.smoothed' | is.null(fs$at) ) stop("Argument must be a result from dlmodeler.filter() or dlmodeler.smooth()")
	if( class(model)!='dlmodeler' | is.null(model$Zt) ) stop("model should be of class 'dlmodeler'")
	if( type[1]!='observation' & type[1]!='state' ) stop("unknown type: ",type[1])
	if( value[1]!='mean' & value[1]!='covariance' & value[1]!='interval' ) stop("unknown value: ",value[1])
	
	if( value[1]=='interval' )
	{
		# extract prediction intervals...
		if( prob<=0 | prob>=1 ) stop("prob should be in the range (0 ; 1)")
		sdfact <- qnorm(.5+prob/2)
		# this could be optimized
		the.mean <- dlmodeler.extract(fs,model,compnames,type,value="mean")
		fun.sd <- function(x) sqrt(apply(x,3,"diag"))
		the.sd <- dlmodeler.extract(fs,model,compnames,type,value="covariance")
		the.sd <- lapply(the.sd,fun.sd)
		ret <- list()
		for(compname in names(the.mean)) {
			ret[[compname]] <- list(
				lower=the.mean[[compname]]-sdfact*the.sd[[compname]],
				mean=the.mean[[compname]],
				upper=the.mean[[compname]]+sdfact*the.sd[[compname]])
		}
		return(ret)
	}
	else
	{
		# which components to extract? if the list is empty, extract
		# everything and give it the model name because we always want
		# this function to return something useful
		if( is.null(compnames) ) components <- model$components else components <- model$components[compnames]
		if( length(components)==0 ) components[model$name] <- matrix(1,NROW(model$Zt),NCOL(model$Zt))
	
		# compute the means and covariances
		# state mean        : E[alpha(t)]   = a'(t)
		# output mean       : E[y(t)]       = Z'(t) %*% a'(t)
		# state covariance  : cov[alpha(t)] = P'(t)
		# output covariance : cov[y(t)]     = Z'(t) %*% P'(t) %*% t(Z'(t)) + H'(t)
		fun <- function(cmp)
		{
			comps <- which(cmp!=0)
			nb.comps <- length(comps)
			n <- NCOL(fs$at) # assuming = dim(fs$Pt)[3]
		
			if( value[1]=='mean' ) {
				# aprim is needed to compute state & output means
				aprim <- matrix(fs$at[comps,],nb.comps,n)
				if( type[1]=='state') return(aprim)
			} else {
				# Pprim is needed to compute state & output covariances
				Pprim <- array(fs$Pt[comps,comps,],dim=c(nb.comps,nb.comps,n))
				if( type[1]=='state') return(Pprim)
			}
			
			# the remainder of this function handles computations for outputs
			# Zprim is needed to compute both output mean & covariance
			d <- NROW(model$Zt)
			if( length(dim(model$Zt))==3 ) {
				Zprim <- array(model$Zt[,comps,],dim=c(d,nb.comps,n))
			} else {
				Zprim <- matrix(model$Zt[,comps],d,nb.comps)
			}
			if( value[1]=='mean' ) {
				# slight modification to aprim
				aprim <- array(aprim,dim=c(nb.comps,1,NCOL(fs$at)))
				ret <- dlmodeler.timevar.fun(Zprim,aprim,match.fun("%*%"))
				return( matrix(ret,d,n) )
			} else {
				fun.zpz <- function(z,p) z %*% p %*% t(z)
				ret <- dlmodeler.timevar.fun(Zprim,Pprim,fun.zpz)
				ret <- dlmodeler.timevar.fun(ret,model$Ht,match.fun("+"))
				return( ret )
			}
		}
		return(lapply(components,fun)) # lapply is magic!
	}
}



dlmodeler.timevar.fun <-
function(x, y, fun)
{
	ndim.x <- length(dim(x))
	ndim.y <- length(dim(y))
	if( ndim.x==1 & ndim.y==1 ) return(fun(x,y)) # ok easy one
	if( ndim.x==2 & ndim.y==2 ) return(fun(x,y)) # also trivial
	if( ndim.x==3 & ndim.y==2 ) {
		# a bit trickier: only one parameter is time-varying
		# deducing the resulting size from the first element
		N <- dim(x)[3]
		# for some strange reason, need to coerce the sub matrix
		# into a matrix...
		guess <- fun(matrix(x[,,1],NROW(x),NCOL(x)),y)
		ret <- array(dim=c(NROW(guess),NCOL(guess),N))
		for( i in 1:N ) ret[,,i] <- fun(matrix(x[,,i],NROW(x),NCOL(x)),y)
		return(ret)
	}
	if( ndim.x==2 & ndim.y==3 ) {
		# same as above
		N <- dim(y)[3]
		guess <- fun(x,matrix(y[,,1],NROW(y),NCOL(y)))
		ret <- array(dim=c(NROW(guess),NCOL(guess),N))
		for( i in 1:N ) ret[,,i] <- fun(x,matrix(y[,,i],NROW(y),NCOL(y)))
		return(ret)
	}
	if( ndim.x==3 & ndim.y==3 ) {
		# same as above, except the two parameters are now time varying as well
		N <- dim(x)[3] # assuming = dim(y)[3]
		guess <- fun(matrix(x[,,1],NROW(x),NCOL(x)),matrix(y[,,1],NROW(y),NCOL(y)))
		ret <- array(dim=c(NROW(guess),NCOL(guess),N))
		for( i in 1:N ) ret[,,i] <- fun(matrix(x[,,i],NROW(x),NCOL(x)),matrix(y[,,i],NROW(y),NCOL(y)))
		return(ret)
	}
	stop("Don't know how to handle dimensions: ",ndim.x,"x",ndim.y)
}



dlmodeler.bdiag <-
function(x, y)
{
	# here is graphically what this function does:
	rbind( cbind(             x             , matrix(0,NROW(x),NCOL(y)) ),
	       cbind( matrix(0,NROW(y),NCOL(x)) ,            y              )
	)
}



dlmodeler.add <-
function(mod1, mod2, name=NULL)
{
	if( is.null(mod1) ) return(mod2)
	if( is.null(mod2) ) return(mod1)
	if( is.null(name) ) name <- paste(mod1$name,mod2$name,sep='+')
	if( class(mod1)!='dlmodeler' ) stop("mod1 should be of class 'dlmodeler'")
	if( class(mod2)!='dlmodeler' ) stop("mod2 should be of class 'dlmodeler'")
	
	# concatenate the state vectors
	a0 <- rbind(mod1$a0,mod2$a0)
	P0 <- dlmodeler.bdiag(mod1$P0,mod2$P0)
	P0inf <- dlmodeler.bdiag(mod1$P0inf,mod2$P0inf)
	Tt <- dlmodeler.timevar.fun(mod1$Tt,mod2$Tt,dlmodeler.bdiag)
	Rt <- dlmodeler.timevar.fun(mod1$Rt,mod2$Rt,dlmodeler.bdiag)
	Qt <- dlmodeler.timevar.fun(mod1$Qt,mod2$Qt,dlmodeler.bdiag)
	
	# add outputs together
	Zt <- dlmodeler.timevar.fun(mod1$Zt,mod2$Zt,cbind)
	myaddition <- function(x,y) (x+y)
	Ht <- dlmodeler.timevar.fun(mod1$Ht,mod2$Ht,myaddition)
	
	# keep track of components
	comp1 <- lapply(mod1$components,function(x) rbind(x,matrix(0,NROW(mod2$Tt),1)) )
	comp2 <- lapply(mod2$components,function(x) rbind(matrix(0,NROW(mod1$Tt),1),x) )
	components <- c(comp1,comp2)
	
	return(dlmodeler.build(a0=a0,P0=P0,P0inf=P0inf,Tt=Tt,Rt=Rt,Qt=Qt,Zt=Zt,Ht=Ht,name=name,components=components))
}



dlmodeler.bind <-
function(mod1, mod2, name=NULL)
{
	if( is.null(mod1) ) return(mod2)
	if( is.null(mod2) ) return(mod1)
	if( is.null(name) ) name <- paste(mod1$name,mod2$name,sep='&')
	if( class(mod1)!='dlmodeler' ) stop("mod1 should be of class 'dlmodeler'")
	if( class(mod2)!='dlmodeler' ) stop("mod2 should be of class 'dlmodeler'")
	
	# concatenate the state vectors
	a0 <- rbind(mod1$a0,mod2$a0)
	P0 <- dlmodeler.bdiag(mod1$P0,mod2$P0)
	P0inf <- dlmodeler.bdiag(mod1$P0inf,mod2$P0inf)
	Tt <- dlmodeler.timevar.fun(mod1$Tt,mod2$Tt,dlmodeler.bdiag)
	Rt <- dlmodeler.timevar.fun(mod1$Rt,mod2$Rt,dlmodeler.bdiag)
	Qt <- dlmodeler.timevar.fun(mod1$Qt,mod2$Qt,dlmodeler.bdiag)
	
	# concatenate the outputs together
	Zt <- dlmodeler.timevar.fun(mod1$Zt,mod2$Zt,dlmodeler.bdiag)
	Ht <- dlmodeler.timevar.fun(mod1$Ht,mod2$Ht,dlmodeler.bdiag)
	
	# keep track of components
	comp1 <- lapply(mod1$components,function(x) rbind(x,matrix(0,NROW(mod2$Tt),1)) )
	comp2 <- lapply(mod2$components,function(x) rbind(matrix(0,NROW(mod1$Tt),1),x) )
	components <- c(comp1,comp2)
	
	return(dlmodeler.build(a0=a0,P0=P0,P0inf=P0inf,Tt=Tt,Rt=Rt,Qt=Qt,Zt=Zt,Ht=Ht,name=name,components=components))
}



 #####################
# modelling functions #
 #####################

dlmodeler.build.polynomial <-
function(ord, sigmaH=1, sigmaQ=0, name='polynomial')
{
	if( ord<0 ) stop("Order must be >= 0")
	m <- ord+1
	if( length(sigmaQ)!=1 & length(sigmaQ)!=m ) stop("SigmaQ has wrong dimension: should be of size ",m)
	d <- 1
	
	a0 <- matrix(0,m,1)
	P0 <- diag(0,m,m)
	P0inf <- diag(m)
	
	Tt <- diag(1,m,m)
	if( m>1 ) for( i in 1:(m-1) ) Tt[i,i+1] <- 1
	Rt <- diag(1,m,m)
	Qt <- diag(sigmaQ^2,m)
	
	Zt <- matrix(c(1,rep(0,m-1)),d,m)
	Ht <- matrix(sigmaH^2,d,d)
	
	return(dlmodeler.build(a0=a0,P0=P0,P0inf=P0inf,Tt=Tt,Rt=Rt,Qt=Qt,Zt=Zt,Ht=Ht,name=name))
}



dlmodeler.build.dseasonal <-
function(ord, sigmaH=1, sigmaQ=0, name='dseasonal')
{
	if( ord<2 ) stop("Order must be >= 2")
	m <- ord-1
	d <- 1
	
	a0 <- matrix(0,m,1)
	P0 <- diag(0,m,m)
	P0inf <- diag(m)
	
	Tt <- matrix(0,m,m)
	if( m>1 ) for( i in 1:(m-1) ) {
		Tt[1,i] <- -1
		Tt[i+1,i] <- 1
	}
	Tt[1,m] <- -1
	Rt <- matrix(0,m,1)
	Rt[1,1] <- 1
	Qt <- matrix(sigmaQ^2,1,1)
	
	Zt <- matrix(c(1,rep(0,m-1)),d,m)
	Ht <- matrix(sigmaH^2,d,d)
	
	return(dlmodeler.build(a0=a0,P0=P0,P0inf=P0inf,Tt=Tt,Rt=Rt,Qt=Qt,Zt=Zt,Ht=Ht,name=name))
}



dlmodeler.build.tseasonal <-
function(per, ord=NULL, sigmaH=1, sigmaQ=0, name='tseasonal')
{
	if( per<=0 ) stop("Period must be > 0")
	if( ((per%%1)!=0) & is.null(ord) ) stop("Order of the trigonometric decomposition must be specified")
	if( is.null(ord) ) m <- per-1 else m <- 2*ord
	if( (m%%1) != 0 ) stop("Order of the trigonopetric decomposition must be an integer value")
	d <- 1
	
	a0 <- matrix(0,m,1)
	P0 <- diag(0,m,m)
	P0inf <- diag(m)
	
	Tt <- diag(-1,m,m)
	if( per==4 ) {
		# bof bof...
		M <- matrix( c(0, -1, 1, 0), 2, 2)
	} else {
		f <- 2*base::pi/per
		M <- matrix( c(cos(f),-sin(f),sin(f),cos(f)), 2, 2 )
	}
	N <- M
	if( m>1) for( i in 1:(m/2) ) {
		Tt[(2*i-1):(2*i), (2*i-1):(2*i)] <- N
		N <- M %*% N
	}
	Rt <- diag(1,m,m)
	Qt <- sigmaQ^2*diag(1,m,m)
	
	Zt <- matrix(rep(c(1,0),m),d,m)
	Ht <- matrix(sigmaH^2,d,d)
	
	return(dlmodeler.build(a0=a0,P0=P0,P0inf=P0inf,Tt=Tt,Rt=Rt,Qt=Qt,Zt=Zt,Ht=Ht,name=name))
}



dlmodeler.build.structural  <-
function(pol.order=NULL, dseas.order=NULL, tseas.order=NULL, tseas.period=NULL,
		sigmaH=1, pol.sigmaQ=0, dseas.sigmaQ=0, tseas.sigmaQ=0, name='structural')
{
	if( !is.null(pol.order) ) {
		if( pol.order==0 ) { pol.name <- 'level'
		} else if( pol.order==1 ) { pol.name <- 'level+trend'
		} else pol.name <- 'polynomial'
		mdl1 <- dlmodeler.build.polynomial(pol.order,sigmaH,pol.sigmaQ,name=pol.name)
	} else mdl1 <- NULL
	if( !is.null(dseas.order) ) mdl2 <- dlmodeler.build.dseasonal(dseas.order,0,dseas.sigmaQ,name='seasonal') else mdl2 <- NULL
	if( !is.null(tseas.order) ) mdl3 <- dlmodeler.build.tseasonal(tseas.period,tseas.order,0,tseas.sigmaQ,name='cycle') else mdl3 <- NULL
	return(dlmodeler.add(mdl1,dlmodeler.add(mdl2,mdl3),name=name))
}



dlmodeler.build.arima <-
function(ar=c(), ma=c(), d=0, sigmaH=1, sigmaQ=0, name='arima')
{
	if( d>0 ) stop("case where d>0 is not implemented yet") # TODO
	if( length(ar)==0 & length(ma)==0 ) stop("ar and/or ma terms are missing")
	r <- max(length(ar),length(ma)+1)
	ar <- c(ar,rep(0,r-length(ar)))
	ma <- c(1,ma,rep(0,r-length(ma)-1))
	
	a0 <- matrix(0,r,1)
	P0 <- diag(0,r,r)
	P0inf <- diag(r)
	
	Tt <- cbind(matrix(ar,r,1), diag(1,r,r-1))
	Rt <- matrix(ma,r,1)
	Qt <- matrix(sigmaQ^2,1,1)
	
	Zt <- matrix(c(1,rep(0,r-1)),1,r)
	Ht <- matrix(sigmaH^2,1,1)
	
	return(dlmodeler.build(a0=a0,P0=P0,P0inf=P0inf,Tt=Tt,Rt=Rt,Qt=Qt,Zt=Zt,Ht=Ht,name=name))
}



dlmodeler.build.regression <-
function(covariates, sigmaH=1, sigmaQ=0, intercept=FALSE, name='regression')
{
	# covariates must be in hoirontal format (1 row per covariate, as y is formatted)
	if( !is.matrix(covariates) ) covariates <- matrix(covariates,nrow=1)
	n <- NCOL(covariates)
	m <- NROW(covariates)+intercept
	if( intercept ) covariates <- rbind(matrix(1,1,n),covariates)
	if( length(sigmaQ)!=1 & length(sigmaQ)!=m ) stop("SigmaQ has wrong dimension: should be of size ",m)
	d <- 1
	a0 <- matrix(0,m,1)
	P0 <- diag(0,m,m)
	P0inf <- diag(m)
	Tt <- diag(1,m,m)
	Rt <- diag(1,m,m)
	Qt <- diag(sigmaQ^2,m,m)
	Zt <- array(dim=c(d,m,n))
	for( i in 1:n ) Zt[,,i] <- t(covariates[,i])
	Ht <- matrix(sigmaH^2,d,d)
	return(dlmodeler.build(a0=a0,P0=P0,P0inf=P0inf,Tt=Tt,Rt=Rt,Qt=Qt,Zt=Zt,Ht=Ht,name=name))
}



 #########################
# filtering and smoothing #
 #########################

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



dlmodeler.filter.FKF <-
function(yt, model, raw.result=FALSE, logLik=TRUE, filter=TRUE)
{
	if(!require('FKF')) stop("required package could not be found: FKF")
	rqr.fun <- function(r,q) r %*% q %*% t(r)
	RQt <- dlmodeler.timevar.fun(model$Rt, model$Qt, rqr.fun)
	m <- NROW(model$Tt)
	d <- NROW(model$Zt)
	res <- FKF::fkf(
			a0=as.vector(model$a0),
			P0=model$P0+1e7*model$P0inf,
			dt=matrix(0,m,1),ct=matrix(0,d,1),
			Tt=model$Tt,
			Zt=model$Zt,
			HHt=RQt,
			GGt=model$Ht,
			yt)
	if( raw.result ) raw.res <- res else raw.res <- NA
	if( length(dim(model$Zt))==2 ) {
		# 1-step ahead prediction when observation matrix is not time-varying
		f <- model$Zt %*% res$at
	} else {
		# 1-step ahead prediction when observation matrix is time-varying
		f <- matrix(NA,NROW(model$Zt),NCOL(yt))
		for( i in 1:NCOL(yt) ) f[,i] <- model$Zt[,,i] %*% res$at[,i]
	}
	return(list(backend='FKF',
					f=f,
					at=res$at,
					Pt=(res$Pt),
					logLik=res$logLik,
					raw.result=raw.res))
}



dlmodeler.filter.dlm <-
function(yt, model, raw.result=FALSE, logLik=FALSE, filter=TRUE)
{
	if(!require('dlm')) stop("required package could not be found: dlm")
	Ttvar <- length(dim(model$Tt))==3
	Rtvar <- length(dim(model$Rt))==3
	Qtvar <- length(dim(model$Qt))==3
	Ztvar <- length(dim(model$Zt))==3
	Htvar <- length(dim(model$Ht))==3
	if( Ttvar|Rtvar|Qtvar|Htvar ) stop("This kind of time varying model is currently unsupported") # TODO handle more time varying cases
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
	return(list(backend='dlm',
					f=res$f,
					at=t(res$a),
					Pt=Pt,
					logLik=logLik,
					raw.result=raw.res))
}



dlmodeler.smooth.dlm <-
function(filt, raw.result=FALSE)
{
	if(!require('dlm')) stop("required package could not be found: dlm")
	res <- dlmSmooth.dlmFiltered(filt$raw.result)
	if( raw.result ) raw.res <- res else raw.res <- NA
	res.Pt <- dlm::dlmSvd2var(res$U.S,res$D.S)
	Pt <- array(NA,dim=c(NCOL(res$D.S),NCOL(res$D.S),length(res.Pt)))
	for( i in 1:length(res.Pt) ) Pt[,,i] <- res.Pt[[i]]
	return(list(backend='dlm',
				at=t(res$s),
				Pt=Pt,
				raw.result=raw.res))
}



dlmodeler.filter <-
function(yt, model, backend=c('KFAS','FKF','dlm'), smooth=FALSE, raw.result=FALSE, logLik=FALSE, filter=TRUE)
{
	if( !is.matrix(yt) ) yt <- matrix(yt,nrow=1)
	if( class(model)!='dlmodeler' ) stop("model should be of class 'dlmodeler'")
	if( backend[1]=='KFAS') { f <- dlmodeler.filter.KFAS(yt,model,raw.result|smooth,logLik,filter)
	} else if( backend[1]=='FKF') { f <- dlmodeler.filter.FKF(yt,model,raw.result|smooth,logLik,filter)
	} else if( backend[1]=='dlm') { f <- dlmodeler.filter.dlm(yt,model,raw.result|smooth,logLik,filter)
	} else stop("Unknown backend: ",backend[1])
	class(f) <- 'dlmodeler.filtered'
	if( smooth ) {
		f$smooth <- dlmodeler.smooth(f,raw.result)
	}
	return(f)
}



dlmodeler.smooth <-
function(filt, raw.result=FALSE)
{
	if( class(filt)!='dlmodeler.filtered' ) stop("argument should be of class 'dlmodeler.filtered'")
	if( is.null(filt$raw.result) ) stop("Raw results are not available; use raw.result=TRUE when calling dlmodeler.filter()")
	backend <- filt$backend
	if( backend[1]=='KFAS') { s <- dlmodeler.smooth.KFAS(filt)
	} else if( backend[1]=='FKF') { stop("FKF does not support smoothing (as of v0.1.1)")
	} else if( backend[1]=='dlm') { s <- dlmodeler.smooth.dlm(filt)
	} else stop("Unknown backend: ",backend[1])
	class(s) <- 'dlmodeler.smoothed'
	return(s)
}



 ###################################
# forecasting and fitting functions #
 ###################################

AIC.dlmodeler.fit <-
function(object, ..., k=2)
{
	if( class(object)!='dlmodeler.fit' ) stop("argument should be of class 'dlmodeler.fit'")
	# AIC pp. 152 of the Durbin & Koopman book
	ll <- object$logLik
	npar <- length(object$par) + sum(object$model$P0inf)
	return( -2*ll + k*npar )
}



dlmodeler.forecast <-
function(yt, model, ahead=1, iters=1, step=1, start=1, prob=.90, backend=c('KFAS','FKF','dlm'))
{
	if( !is.matrix(yt) ) yt <- matrix(yt,nrow=1)
	if(NROW(yt)>1) stop("Multivariate case not implemented yet") # TODO
	if( class(model)!='dlmodeler' ) stop("model should be of class '")
	if( ahead<1 ) stop("ahead should be >= 1")
	if( iters<1 ) stop("iters should be >= 1")
	if( step<1 ) stop("step should be >= 1")
	if( start<1 ) stop("start should be >= 1")
	if( prob<=0 | prob>=1 ) stop("prob should be in the range (0 ; 1)")
	
	sdfact <- qnorm(.5+prob/2)
	mymodel <- model
	pos <- start
	obs <- matrix(yt[1,1:pos],1,pos)
	if(ahead>1) nas <- matrix(NA,1,ahead-1) else nas <- NULL
	simus <- data.frame()
	
	# TODO: maybe we want to be able to forecast for states
	for(i in 1:iters)
	{
		obs <- cbind(obs,nas)
		modfilter <- dlmodeler.filter(obs,mymodel,backend)
		modfilter.std <- dlmodeler.extract(modfilter,mymodel,mymodel$name,type="observation",value="covariance")[[mymodel$name]]
		modfilter.std <- matrix(sqrt(modfilter.std),nrow=1)
		
		# FIXME f is already a 1-step ahead forecast
		sim <- data.frame(
				index=(pos+1):(pos+ahead),
				distance=1:ahead,
				lower=t(modfilter$f)[(pos+1):(pos+ahead)]-sdfact*t(modfilter.std)[(pos+1):(pos+ahead)],
				yhat=t(modfilter$f)[(pos+1):(pos+ahead)],
				upper=t(modfilter$f)[(pos+1):(pos+ahead)]+sdfact*t(modfilter.std)[(pos+1):(pos+ahead)],
				y=yt[(pos+1):(pos+ahead)])
		simus <- rbind(simus,sim)
		
#		TODO: restart filtering from last step to speed up the simulation
#		last <- nrow(modfilter$m)
#		if(NCOL(modfilter$m)>0) {
#			mymodel$m0 <- modfilter$m[last,]
#			mymodel$C0 <- with(modfilter,dlmSvd2var(U.C[[last]],D.C[last,]))
#		} else {
#			mymodel$m0 <- modfilter$m[last]
#			mymodel$C0 <- with(modfilter,dlmSvd2var(U.C,D.C[last]))
#		}
#		pos <- pos+step
#		obs <- yt[(pos-step+1):(pos)]
		
		pos <- pos+step
		obs <- matrix(yt[1:pos],1,pos)
	}
	
	return(simus)
}



dlmodeler.fit.MLE <-
function(yt, build.fun, par, backend=c('KFAS','FKF','dlm'),
		method="L-BFGS-B", verbose=FALSE,
		filter=TRUE, smooth=FALSE, raw.result=FALSE, ...)
{
	fit.fun <- function(p, ...) {
		f <- dlmodeler.filter(yt, build.fun(p), backend=backend, logLik=TRUE, filter=FALSE)
		if( verbose ) cat(p,':',-f$logLik,'\n')
		return(-f$logLik)
	}
	opt <- optim(par, fit.fun, ..., method=method)
	if( verbose ) cat(opt$message,":",opt$convergence,"\n")
	opt$model <- build.fun(opt$par)
	if(filter) opt$filtered <- dlmodeler.filter(yt, opt$model, backend=backend, smooth=smooth, raw.result=raw.result, logLik=TRUE)
	opt$logLik <- -opt$value
	opt$par0 <- par
	class(opt) <- "dlmodeler.fit"
	return(opt)
}



dlmodeler.fit.MSE <-
function(yt, build.fun, par, ahead, iters=NCOL(yt)-ahead, step=1, start=1,
		backend=c('KFAS','FKF','dlm'), method="L-BFGS-B", verbose=FALSE,
		filter=TRUE, smooth=FALSE, raw.result=FALSE, ...)
{
	fit.fun <- function(p, ...) {
		f <- dlmodeler.forecast(yt, build.fun(p), backend=backend,
				ahead=ahead, iters=iters, step=step, start=start)
		mse <- mean((f$yhat-f$y)^2,na.rm=TRUE)
		if( verbose ) cat(p,':',mse,'\n')
		return(mse)
	}
	opt <- optim(par, fit.fun, ..., method=method)
	if( verbose ) cat(opt$message,":",opt$convergence,"\n")
	opt$model <- build.fun(opt$par)
	if(filter) opt$filtered <- dlmodeler.filter(yt, opt$model, backend=backend, smooth=smooth, raw.result=raw.result, logLik=TRUE)
	opt$logLik <- opt$filtered$logLik
	opt$par0 <- par
	class(opt) <- "dlmodeler.fit"
	return(opt)
}	



dlmodeler.fit.MAD <-
function(yt, build.fun, par, ahead, iters=NCOL(yt)-ahead, step=1, start=1,
		backend=c('KFAS','FKF','dlm'), method="L-BFGS-B", verbose=FALSE,
		filter=TRUE, smooth=FALSE, raw.result=FALSE, ...)
{
	fit.fun <- function(p, ...) {
		f <- dlmodeler.forecast(yt, build.fun(p), backend=backend,
				ahead=ahead, iters=iters, step=step, start=start)
		mse <- mean(abs(f$yhat-f$y),na.rm=TRUE)
		if( verbose ) cat(p,':',mse,'\n')
		return(mse)
	}
	opt <- optim(par, fit.fun, ..., method=method)
	if( verbose ) cat(opt$message,":",opt$convergence,"\n")
	opt$model <- build.fun(opt$par)
	if(filter) opt$filtered <- dlmodeler.filter(yt, opt$model, backend=backend, smooth=smooth, raw.result=raw.result, logLik=TRUE)
	opt$logLik <- opt$filtered$logLik
	opt$par0 <- par
	class(opt) <- "dlmodeler.fit"
	return(opt)
}	



dlmodeler.fit.sigma <-
function(yt, build.fun, par, backend=c('KFAS','FKF','dlm'), method="L-BFGS-B",
		verbose=FALSE, filter=TRUE, smooth=FALSE, raw.result=FALSE, ...)
{
	fit.fun <- function(p, ...) {
		f <- dlmodeler.filter(yt, build.fun(p), backend=backend)
		sigma <- sqrt(var(f$f[1,1:NCOL(yt)] - yt[1,]))
		if( verbose ) cat(p,':',sigma,'\n')
		return(sigma)
	}
	opt <- optim(par, fit.fun, ..., method=method)
	if( verbose ) cat(opt$message,":",opt$convergence,"\n")
	opt$model <- build.fun(opt$par)
	if(filter) opt$filtered <- dlmodeler.filter(yt, opt$model, backend=backend, smooth=smooth, raw.result=raw.result, logLik=TRUE)
	opt$logLik <- opt$filtered$logLik
	opt$par0 <- par
	class(opt) <- "dlmodeler.fit"
	return(opt)
}





