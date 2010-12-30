
require(dlmodeler)

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

# print the model
print(mod)

# check if the model is valid
dlmodeler.check(mod)[1]==1
dlmodeler.check(mod)[2]==2
dlmodeler.check(mod)[3]==2
dlmodeler.check(mod)[4]==1
dlmodeler.check(mod)[5]==0
is.na(dlmodeler.check(mod)[6])
is.na(dlmodeler.check(mod)[7])
is.na(dlmodeler.check(mod)[8])
is.na(dlmodeler.check(mod)[9])
is.na(dlmodeler.check(mod)[10])


# an empty DLM with 4 state variables (3 of which are stocastic)
# and bi-variate observations
mod <- dlmodeler.build(dimensions=c(4,3,2))

# check if the model is valid
# check if the model is valid
dlmodeler.check(mod)[1]==1
dlmodeler.check(mod)[2]==4
dlmodeler.check(mod)[3]==3
dlmodeler.check(mod)[4]==2
dlmodeler.check(mod)[5]==0
is.na(dlmodeler.check(mod)[6])
is.na(dlmodeler.check(mod)[7])
is.na(dlmodeler.check(mod)[8])
is.na(dlmodeler.check(mod)[9])
is.na(dlmodeler.check(mod)[10])

# print the model
print(mod)

# this example is fairly complete, test it with various backends
# generate some data
N <- 365*5
t <- c(1:N,rep(NA,365))
a <- rnorm(N+365,0,.5)
y <- pi + cos(2*pi*t/365.25) + .25*sin(2*pi*t/365.25*3) +
		exp(1)*a + rnorm(N+365,0,.5)

# build a model for this data
m1 <- dlmodeler.build.polynomial(0,sigmaH=.5,name='level')
m2 <- dlmodeler.build.dseasonal(7,sigmaH=0,name='week')
m3 <- dlmodeler.build.tseasonal(365.25,3,sigmaH=0,name='year')
m4 <- dlmodeler.build.regression(a,sigmaH=0,name='reg')
m <- dlmodeler.add(m1,
		dlmodeler.add(m2,dlmodeler.add(m3,m4)),
		name='mymodel')

test.backend <- function(backend)
{
	cat(backend,"\n")
	
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
	m.state.mean$level[,N] # should be close to pi
	mean(abs(m.state.mean$week[,N])) # should be close to 0
	m.state.mean$year[1,N] # should be close to 1
	m.state.mean$year[6,N] # should be close to .25
	m.state.mean$reg[,N] # should be close to e
	
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
	
	return(TRUE)
}

backs = list(KFAS='KFAS',FKF='FKF',dlm='dlm')
lapply(backs, test.backend)

