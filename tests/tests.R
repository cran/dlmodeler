
require(dlmodeler)

####################
# tests for 'core' #
####################

# create a DLM by specifying its vectors and matrices
# check if the model is valid
mod <- dlmodeler.build(
		a0 = c(0,0), # initial state: (level, trend)
		P0 = diag(c(0,0)), # initial state variance set to...
		P0inf = diag(2), # ...use exact diffuse initialization
		matrix(c(1,0,1,1),2,2), # state transition matrix
		diag(c(1,1)), # state disturbance selection matrix
		diag(c(.5,.05)), # state disturbance variance matrix
		matrix(c(1,0),1,2), # observation design matrix
		matrix(1,1,1) # observation disturbance variance matrix
)

print(mod)
dlmodeler.check(mod)$status==TRUE
dlmodeler.check(mod)$m==2
dlmodeler.check(mod)$r==2
dlmodeler.check(mod)$d==1
dlmodeler.check(mod)$timevar==FALSE
is.na(dlmodeler.check(mod)$timevar.Tt)
is.na(dlmodeler.check(mod)$timevar.Rt)
is.na(dlmodeler.check(mod)$timevar.Qt)
is.na(dlmodeler.check(mod)$timevar.Zt)
is.na(dlmodeler.check(mod)$timevar.Ht)


# an empty DLM with 4 state variables (3 of which are stochastic)
# and bi-variate observations, check if the model is valid
mod <- dlmodeler.build(dimensions=c(4,3,2))

print(mod)
dlmodeler.check(mod)$status==TRUE
dlmodeler.check(mod)$m==4
dlmodeler.check(mod)$r==3
dlmodeler.check(mod)$d==2
dlmodeler.check(mod)$timevar==FALSE
is.na(dlmodeler.check(mod)$timevar.Tt)
is.na(dlmodeler.check(mod)$timevar.Rt)
is.na(dlmodeler.check(mod)$timevar.Qt)
is.na(dlmodeler.check(mod)$timevar.Zt)
is.na(dlmodeler.check(mod)$timevar.Ht)


# operations on matrices
v1 <- matrix(1:9,nrow=3,ncol=3)
v2 <- array(1:18,dim=c(3,3,2))

m1 <- dlmodeler.timevar.fun(v1,v1,function(x,y) x+y)
sum(abs(m1-v1-v1))==0
m21 <- dlmodeler.timevar.fun(v2,v1,function(x,y) x+y)
sum(abs(m21[,,1]-v1-v1))==0
sum(abs(m21[,,2]-v1-v1-9))==0
m22 <- dlmodeler.timevar.fun(v1,v2,function(x,y) x+y)
sum(abs(m22[,,1]-m21[,,1]))==0
sum(abs(m22[,,2]-m21[,,2]))==0
m3 <- dlmodeler.timevar.fun(v2,v1,function(x,y) x+y)
sum(abs(m3[,,1]-v1-v1))==0
sum(abs(m3[,,2]-v1-matrix(10:18,nrow=3,ncol=3)))==0

mt <- matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,9,9,0,0,0,9,9,9,0,0,0,9,9,9),nrow=6,ncol=6)
md <- dlmodeler.timevar.fun(v1,v2,dlmodeler.bdiag)
sum(abs(md[,,2]-md[,,1]-mt))==0


# this example is fairly complete, covers 'add', 'filter', 'smooth'
# 'extract', 'polynomial', 'dseasonal', 'tseasonal', 'regression'
# test it with various backends
# generate some data
set.seed(19820605)
N <- 365*10
t <- c(1:N,rep(NA,365))
a <- rnorm(N+365,0,.5)
y <- pi + cos(2*pi*t/365.25) + .25*sin(2*pi*t/365.25*3) +
		exp(1)*a + rnorm(N+365,0,.5)

# build a model for this data
m1 <- dlmodeler.build.polynomial(0,sigmaH=.5,name='level')
m2 <- dlmodeler.build.dseasonal(7,sigmaH=0,name='week')
m3 <- dlmodeler.build.tseasonal(365.25,3,sigmaH=0,name='year')
m4 <- dlmodeler.build.regression(a,sigmaH=0,name='reg')
m <- dlmodeler.add(m1,dlmodeler.add(m2,dlmodeler.add(m3,m4)),
		name='mymodel')

test.backend <- function(backend)
{
	cat(backend,"\n")
	test.ok <- TRUE
	
	system.time(f <- dlmodeler.filter(y, m, raw.result=TRUE, backend=backend))
	
	# extract all the components
	m.state.mean <- dlmodeler.extract(f,m,type="state",value="mean")
	m.state.cov <- dlmodeler.extract(f,m,type="state",value="covariance")
	m.obs.mean <- dlmodeler.extract(f,m,type="observation",value="mean")
	m.obs.cov <- dlmodeler.extract(f,m,type="observation",value="covariance")
	m.obs.int <- dlmodeler.extract(f,m,type="observation",value="interval",prob=.01)
	
	par(mfrow=c(2,1))
	
	# show the one step ahead forecasts & 99\% prediction intervals
	plot(y,xlim=c(N-10,N+30))
	lines(m.obs.int$mymodel$upper[1,],col='light grey')
	lines(m.obs.int$mymodel$lower[1,],col='light grey')
	lines(m.obs.int$mymodel$mean[1,],col=2)
	
	# see to which values the filter has converged:
	test.ok <- test.ok & abs(m.state.mean$level[,N]-pi)/pi < .05 # should be close to pi
	test.ok <- test.ok & abs(mean(abs(m.state.mean$week[,N]))) < .05 # should be close to 0
	test.ok <- test.ok & abs(m.state.mean$year[1,N]-1) < .05 # should be close to 1
	test.ok <- test.ok & abs(m.state.mean$year[6,N]-.25) < .05 # should be close to .25
	test.ok <- test.ok & abs(m.state.mean$reg[,N]-exp(1))/exp(1) < .05 # should be close to e
	
	# show the filtered level+year components
	plot(m.obs.mean$level[1,]+m.obs.mean$year[1,],
			type='l',ylim=c(pi-2,pi+2),col='light green',
			ylab="smoothed & filtered level+year")
	
	if(backend!='FKF') {
		system.time(s <- dlmodeler.smooth(f))
	
		# show the smoothed level+year components
		s.obs.mean <- dlmodeler.extract(s,m,type="observation",value="mean")
		lines(s.obs.mean$level[1,]+s.obs.mean$year[1,],type='l',
			ylim=c(pi-2,pi+2),col='dark green')
	}
	
	return(test.ok)
}

test.backend('KFAS')
test.backend('FKF')
test.backend('dlm')



##################################
# tests for 'dlm', 'FKF', 'KFAS' #
##################################

yt <- matrix(c(Nile,rep(NA,12)),nrow=1)
mod <- dlmodeler.build.polynomial(1,sigmaH=1,sigmaQ=c(exp(-5),exp(-8)),name='mymodel')
mod$a0 <- c(10,40)
mod$P0 <- mod$P0inf*1e2
mod$P0inf <- mod$P0inf*0

ref <- dlmodeler.filter(yt,mod,backend='KFAS')
test.backend <- function(backend)
{
	mod$a0 <- c(10,40)
	filt <- dlmodeler.filter(c(yt),mod,backend=backend)
	mean(abs(filt$f[1,5:100]-ref$f[1,5:100]))
}

test.backend('KFAS')/mean(yt,na.rm=TRUE) < .05
test.backend('FKF')/mean(yt,na.rm=TRUE) < .05
test.backend('dlm')/mean(yt,na.rm=TRUE) < .05


#############################
# advanced tests for 'core' #
#############################

# forecast simulation with a constant model (no diffuse init)
yt <- matrix(c(Nile,rep(NA,12)),nrow=1)
mod <- dlmodeler.build.polynomial(1,sigmaH=1,sigmaQ=c(exp(-5),exp(-8)),name='mymodel')
mod$a0 <- c(0,40)
mod$P0 <- mod$P0inf*1e3
mod$P0inf <- mod$P0inf*0

test.backend <- function(backend)
{
	nb.ahead <- 5
	nb.iters <- 10
	nb.start <- 10
	fcst.debug <- dlmodeler.forecast(yt, mod, nb.ahead, nb.iters, start=nb.start, backend=backend, debug=TRUE)
	fcst.fast <- dlmodeler.forecast(yt, mod, nb.ahead, nb.iters, start=nb.start, backend=backend, debug=FALSE)
	ctl <- dlmodeler.filter(yt,mod,backend='KFAS')$f[1,nb.start:(nb.iters+nb.start-1)]
	ctl.debug <- subset(fcst.debug,distance==1)$yhat
	ctl.fast <- subset(fcst.fast,distance==1)$yhat
	d0.fast <- sum(abs(ctl-ctl.debug))
	d0.debug <- sum(abs(ctl-ctl.fast))

	nb.ahead <- 5
	nb.iters <- 10
	nb.start <- 50
	fcst.debug <- dlmodeler.forecast(yt, mod, nb.ahead, nb.iters, start=nb.start, backend=backend, debug=TRUE)
	fcst.fast <- dlmodeler.forecast(yt, mod, nb.ahead, nb.iters, start=nb.start, backend=backend, debug=FALSE)
	ctl <- dlmodeler.filter(yt,mod,backend='KFAS')$f[1,nb.start:(nb.iters+nb.start-1)]
	ctl.debug <- subset(fcst.debug,distance==1)$yhat
	ctl.fast <- subset(fcst.fast,distance==1)$yhat
	d1.fast <- sum(abs(ctl-ctl.debug))
	d1.debug <- sum(abs(ctl-ctl.fast))
	
	return(d0.debug+d0.fast+d1.debug+d1.fast)
}

test.backend('KFAS')/mean(yt,na.rm=TRUE) < 0.05
test.backend('FKF')/mean(yt,na.rm=TRUE) < 0.05
test.backend('dlm')/mean(yt,na.rm=TRUE) < 5

# forecast simulation with a constant model (with diffuse)
yt <- matrix(c(Nile,rep(NA,12)),nrow=1)
mod <- dlmodeler.build.polynomial(1,sigmaH=1,sigmaQ=c(exp(-5),exp(-8)),name='mymodel')
mod$a0 <- c(0,40)

test.backend('KFAS')/mean(yt,na.rm=TRUE) < 0.05
test.backend('FKF')/mean(yt,na.rm=TRUE) < 0.05
test.backend('dlm')/mean(yt,na.rm=TRUE) < 5

# TODO: tests for 'bind'



###################
# tests for 'fit' #
###################

y <- matrix(Nile,nrow=1)
build.fun <- function(p) {
	varH <- exp(p[1])
	varQ <- exp(p[2])
	dlmodeler.build.polynomial(0,sqrt(varH),sqrt(varQ),name='p32')
}

# fit the model by maximum likelihood estimation
# compare the fitted parameters with those reported by the Durbin & Koopman
fit <- dlmodeler.fit.MLE(y, build.fun, c(0,0), verbose=FALSE)
abs(fit$model$Ht-15099)/15099 < .05
abs(fit$model$Qt-1469.1)/1469.1 < .05

# test build.unknowns
test.model <- build.fun(c(NA,NA))
test.model <-dlmodeler.build.unknowns(test.model)
test.model
length(test.model$unknowns)==2

# test build.function
bfun <- dlmodeler.build.function(test.model)
my.model <- bfun(fit$par)
my.model
length(my.model$unknowns)==0
abs(my.model$Ht-15099)/15099 < .05
abs(my.model$Qt-1469.1)/1469.1 < .05

# test automatic fitting with NA
fit <- dlmodeler.fit(yt,build.fun(c(NA,NA)),method="MLE")
(fit$model$Ht-15099)/15099 < .05
(fit$model$Qt-1469.1)/1469.1 < .05

# test AIC & logLik
logLik(fit)
abs(-2*logLik(fit)+2*3-AIC(fit))/AIC(fit) < .05


# test fitting functions
fit.MLE <- dlmodeler.fit(yt,build.fun(c(NA,NA)),method="MLE")
fit.MSE <- dlmodeler.fit(yt,build.fun(c(NA,NA)),method="MSE",ahead=5,start=20)
fit.MAD <- dlmodeler.fit(yt,build.fun(c(NA,NA)),method="MAD",ahead=5,start=20)
fit.MAPE <- dlmodeler.fit(yt,build.fun(c(NA,NA)),method="MAPE",ahead=1,start=20)

fcst.MLE <- dlmodeler.forecast(yt, fit.MLE$model, ahead=4, iters=40, start=50)
fcst.MSE <- dlmodeler.forecast(yt, fit.MSE$model, ahead=4, iters=40, start=50)
fcst.MAD <- dlmodeler.forecast(yt, fit.MAD$model, ahead=4, iters=40, start=50)
fcst.MAPE <- dlmodeler.forecast(yt, fit.MAPE$model, ahead=4, iters=40, start=50)

mean(abs(fcst.MLE$yhat-fcst.MLE$y)) < 100
mean(abs(fcst.MSE$yhat-fcst.MSE$y)) < 100
mean(abs(fcst.MAD$yhat-fcst.MAD$y)) < 100
mean(abs(fcst.MAPE$yhat-fcst.MAPE$y)) < 100


#
# tests for 'models-basic'
#

# TODO: tests for 'structural', 'regression', 'arima'



#
# tests for 'models-advanced'
#

# TODO: tests for 'weekdays'


