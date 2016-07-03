

## tracking down error in conjugacy checking system,
## this model was submitted by: Eduardo Martins <egmartins@gmail.com>
## during the ISEC 2016 NIMBLE workshop
library(nimble)
load('~/Downloads/data_inits.RData')
fec_code <- nimbleCode({
    ## Likelihood
    ## Fecundity
    for (i in 1:ny.fec) {
        y.fec[i] ~ dlnorm(log(z.fec[i]), tau.y.fec) # observation model
        z.fec[i] ~ dnegbin(p.fec[i], size.fec) # sampling model
        p.fec[i] <- size.fec / (size.fec + mu.fec[i])
        log(mu.fec[i]) <- a.fec.raw + b.fec * z.sl[i] + nu.fec.k.raw[idx.k[i]] + nu.fec.l.raw[idx.l[i]]
    }
    a.fec <- a.fec.raw + mean(nu.fec.k.raw[1:nk]) + mean(nu.fec.l.raw[1:nl])
    for (k in 1:nk){
        nu.fec.k.raw[k] ~ dnorm(mu.nu.fec.k, tau.nu.fec.k)
        nu.fec.k[k] <- nu.fec.k.raw[k] - mean(nu.fec.k.raw[1:nk])
    }
    for (l in 1:nl){
        nu.fec.l.raw[l] ~ dnorm(mu.nu.fec.l, tau.nu.fec.l)
        nu.fec.l[l] <- nu.fec.l.raw[l] - mean(nu.fec.l.raw[1:nl])
    }
    ## Length
    for(i in 1:ny.fec) {
        z.sl[i] ~ dnorm(mu.sl[i], tau.sl.i)
        mu.sl[i] <- a.sl.raw + b.sl * z.ocr[i] + nu.sl.k.raw[idx.k[i]] + nu.sl.l.raw[idx.l[i]]
    }
    a.sl <- a.sl.raw + mean(nu.sl.k.raw[1:nk]) + mean(nu.sl.l.raw[1:nl])
    for (k in 1:nk){
        nu.sl.k.raw[k] ~ dnorm(mu.nu.sl.k, tau.nu.sl.k)
        nu.sl.k[k] <- nu.sl.k.raw[k] - mean(nu.sl.k.raw[1:nk])
    }
    for (l in 1:nl){
        nu.sl.l.raw[l] ~ dnorm(mu.nu.sl.l[l], tau.nu.sl.l)
        mu.nu.sl.l[l] <- b.nu.sl.l[1] * sst6.std[l] + b.nu.sl.l[2] * pdo6.std[l] + b.nu.sl.l[3] * nsal.std[l]
        nu.sl.l[l] <- nu.sl.l.raw[l] - mean(nu.sl.l.raw[1:nl])
    }
    ## Years of ocean residence
    for(i in 1:ny.fec){
        z.ocr[i] ~ dbern(p.z.ocr[i])
        logit(p.z.ocr[i]) <- a.ocr.raw + nu.ocr.k.raw[idx.k[i]] + nu.ocr.l.raw[idx.l[i]]
    }
    a.ocr <- a.ocr.raw + mean(nu.ocr.k.raw[1:nk]) + mean(nu.ocr.l.raw[1:nl])
    for (k in 1:nk){
        nu.ocr.k.raw[k] ~ dnorm(mu.nu.ocr.k, tau.nu.ocr.k)
        nu.ocr.k[k] <- nu.ocr.k.raw[k] - mean(nu.ocr.k.raw[1:nk])
    }
    for (l in 1:nl){
        nu.ocr.l.raw[l] ~ dnorm(mu.nu.ocr.l, tau.nu.ocr.l)
        nu.ocr.l[l] <- nu.ocr.l.raw[l] - mean(nu.ocr.l.raw[1:nl])
    }
    ## Priors
    ## Fecundity
    tau.y.fec <- pow(sig.y.fec, -2)
    sig.y.fec ~ dunif(0, 100)
    num.size.fec ~ dnorm(0, 0.0016)
    den.size.fec ~ dnorm(0, 1)
    size.fec <- abs(num.size.fec / den.size.fec)
    mu.nu.fec.k ~ dnorm(0, 0.0001)
    tau.nu.fec.k <- pow(sig.nu.fec.k, -2)
    sig.nu.fec.k ~ dunif(0, 100)
    mu.nu.fec.l ~ dnorm(0, 0.0001)
    tau.nu.fec.l <- pow(sig.nu.fec.l, -2)
    sig.nu.fec.l ~ dunif(0, 100)
    a.fec.raw ~ dnorm(mean.z.fec, tau.z.fec)
    b.fec ~ dnorm(0, 0.0001)
    ## Length
    tau.sl.i <- pow(sig.sl.i, -2)
    sig.sl.i ~ dunif(0, 100)
    mu.nu.sl.k ~ dnorm(0, 0.0001)
    tau.nu.sl.k <- pow(sig.nu.sl.k, -2)
    sig.nu.sl.k ~ dunif(0, 100)
    tau.nu.sl.l <- pow(sig.nu.sl.l, -2)
    sig.nu.sl.l ~ dunif(0, 100)
    a.sl.raw ~ dnorm(0, 0.0001)
    b.sl ~ dnorm(0, 0.0001)
    for (b in 1:nb.nu.sl.l) {
        b.nu.sl.l[b] ~ dnorm(0, 0.0001)
    }
    ## Years of ocean residence
    mu.nu.ocr.k ~ dnorm(0, 0.0001)     
    tau.nu.ocr.k <- pow(sig.nu.ocr.k, -2)
    sig.nu.ocr.k ~ dunif(0, 100)
    mu.nu.ocr.l ~ dnorm(0, 0.0001)
    tau.nu.ocr.l <- pow(sig.nu.ocr.l, -2)
    sig.nu.ocr.l ~ dunif(0, 100)
    a.ocr.raw ~ dnorm(0, 0.0001)
    ## Residuals and metrics for model assessment
    ## Fecundity
    for(i in 1:ny.fec){
        res.y.fec[i] <- y.fec[i] - z.fec[i]
        new.y.fec[i] ~ dlnorm(log(z.fec[i]), tau.y.fec)
    }
    for(i in 1:ny.fec){
        pres.z.fec[i] <- (z.fec[i] - mu.fec[i])/sqrt(mu.fec[i] + pow(mu.fec[i], 2) / size.fec)
        new.z.fec[i] ~ dnegbin(p.fec[i], size.fec)
        pres.new.z.fec[i] <- (new.z.fec[i] - mu.fec[i])/sqrt(mu.fec[i] + pow(mu.fec[i], 2) / size.fec)
        d.z.fec[i] <- pow(pres.z.fec[i], 2)
        d.new.z.fec[i] <- pow(pres.new.z.fec[i], 2)
    }
    fit.z.fec <- sum(d.z.fec[1:ny.fec])
    fit.new.z.fec <- sum(d.new.z.fec[1:ny.fec])
    ## Length
    for (i in 1:ny.fec) {
        res.z.sl[i] <- z.sl[i] - mu.sl[i]
        new.z.sl[i] ~ dnorm(mu.sl[i], tau.sl.i)
    }
    for (l in 1:nl){
        res.nu.sl.l[l] <- nu.sl.l[l] - mu.nu.sl.l[l]
        new.nu.sl.l[l] ~ dnorm(mu.nu.sl.l[l], tau.nu.sl.l)
    }
})
fec_mod <- nimbleModel(fec_code, constants = dat, inits = inits())
fec_cmod <- compileNimble(fec_mod)

##options(error = recover)

## next line errors out in conjugacy check:
## Error in x[[1]] : subscript out of bounds
mcmcConf <- configureMCMC(fec_mod)

mcmcConf$printSamplers()
mcmcConf$getMonitors()
fecMCMC <- buildMCMC(mcmcConf)
cfecMCMC <- compileNimble(fecMCMC, project = fec_mod)

cfecMCMC$run(10)

## testing of new setData() functionality:
## taking a character vector of variable names

library(nimble)
Rmodel <- readBUGSmodel('birats2.bug', dir = getBUGSexampleDir('birats'), data = 'birats-data.R', inits = 'birats-inits.R')

Rmodel$getNodeNames(includeData = FALSE)

modelLs <- readBUGSmodel('birats2.bug', dir = getBUGSexampleDir('birats'), data = 'birats-data.R', inits = 'birats-inits.R', returnModelComponents = TRUE)
modelLs$model

Rmodel$beta
Rmodel$tau.c
Rmodel$r

Rmodel$getNodeNames(dataOnly = TRUE)
Rmodel$setData(c('beta', 'tau.c', 'r'))
Rmodel$getNodeNames(dataOnly = TRUE)
Rmodel$resetData()
Rmodel$getNodeNames(dataOnly = TRUE)
Rmodel$setData(c('beta', 'tau.c', 'r'))
Rmodel$getNodeNames(dataOnly = TRUE)
Rmodel$setData(c('Y'))





## testing new MCMC runtime argument: time
library(nimble)

code <- nimbleCode({
    a ~ dnorm(0, 1)
    b ~ dnorm(a, 1)
    c ~ dnorm(b + b^2, 1)
})
constants <- list()
data <- list(b=0, c=0)
inits <- list(a=0)

Rmodel <- nimbleModel(code, constants, data, inits)

spec <- configureMCMC(Rmodel)
spec$printSamplers()
Rmcmc <- buildMCMC(spec)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Cmcmc$run(100000, time = TRUE)
Cmcmc$getTimes()
Cmcmc$run(100000, time = TRUE, reset = FALSE)
Cmcmc$getTimes()



## error from Perry about not finding nimArray
## but I can't reproduce this error

library(nimble)

Rmodel <- readBUGSmodel('birats2.bug', dir = getBUGSexampleDir('birats'), data = 'birats-data.R', inits = 'birats-inits.R')
##spec <- configureMCMC(Rmodel)
##spec$printSamplers()
##spec$getSamplerDefinition(1)
##spec$getSamplerDefinition(4)
##mcmc <- buildMCMC(spec)
mcmc <- buildMCMC(Rmodel)
mcmc$run(1)
##Error in samplerFunctions[[i]]$run() : could not find function "nimArray"

source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))
test_mcmc('birats', model = 'birats2.bug', inits = 'birats-inits.R', data = 'birats-data.R', numItsC = 1000, resampleData = TRUE)




## Soledad's compilation error

library(nimble)
x <- c(2, 2, 3)

## works
test1 <- nimbleFunction ( run = function (indx = double(1) ) {
    if(  length(indx)==2 & indx[1]==indx[2] ) {
        aa <- indx[1] 
    } else {
        newindx <- indx
        aa         <- length(newindx)
    }
    output <- 1 + aa
    returnType(double(0))
    return(output)
})

## doesn't work
test1 <- nimbleFunction(run = function(indx = double(1)) {
    if(  length(indx)==2 & indx[1]==indx[2] ) {
        newindx <- indx[1]
        aa          <- newindx     ## THIS IS THE MAIN CHANGE IN THIS FUNCTION
    } else {
        newindx <- indx
        aa          <- length(newindx)
    }
    output <- 1 + aa
    returnType(double(0))
    return(output)
})

Ctest1 <- compileNimble(test1)


test1(x)
Ctest1(x)





## using a custom distribution to have flexible length data
## for floriane plard population models

library(nimble)


dxxx <- nimbleFunction(
    run = function(x = double(1), mu = double(1), sigma = double(1), length = integer(), log.p = double()) {
        ll <- 0
        for(i in 1:length) {
            ll <- ll + dnorm(x[i], mu[i], sd=sigma[i], log=TRUE)
        }
        returnType(double())
        if(log.p) return(ll) else return(exp(ll))
    }
)

rxxx <- nimbleFunction(
    run = function(n = integer(), mu = double(1), sigma = double(1), length = integer()) {
        print('this should never run')
        x <- numeric(length)
        return(x)
    }
)

registerDistributions(list(
    dxxx = list(
        BUGSdist = 'dxxx(mu, sigma, length)',
        types    = c('value = double(1)', 'mu = double(1)', 'sigma = double(1)', 'length = integer()')
    )
))

code <- nimbleCode({
    y[1:N] ~ dxxx(mu[1:N], sigma[1:N], length)
})
constants <- list(N = 10)
data <- list()
inits <- list(length = 10)

Rmodel <- nimbleModel(code, constants, data, inits)

Rmodel$length
Rmodel$mu
Rmodel$sigma
Rmodel$y
Rmodel$getNodeNames(dataOnly = TRUE)

length <- 10
length
mu <- ((1:10)/2)[1:length]
mu
sigma <- (1+(1:10)/10)[1:length]
sigma
y <- c(1,0,1,3,2,5,4,6,7,5)[1:length]
y
yaug <- c(y, rep(as.numeric(NA), 10-length))
yaug

Rmodel$length <- length
Rmodel$length
Rmodel$mu[1:length] <- mu
Rmodel$mu
Rmodel$sigma[1:length] <- sigma
Rmodel$sigma
Rmodel$resetData()
Rmodel$getNodeNames(dataOnly = TRUE)
Rmodel$setData(list(y=yaug))
Rmodel$getNodeNames(dataOnly = TRUE)
Rmodel$y

sum(dnorm(y, mu, sd=sigma, log=TRUE))
calculate(Rmodel)


debug(exprClasses_setSizes)



## testing MCEM for pump model
library(nimble)

pumpCode <- nimbleCode({ 
  for (i in 1:N){
      theta[i] ~ dgamma(alpha,beta) 
      lambda[i] <- theta[i]*t[i]
      x[i] ~ dpois(lambda[i]) 
  }
  alpha ~ dexp(1.0) 
  beta ~ dgamma(0.1,1.0) 
}) 
pumpConsts <- list(N = 10,
                   t = c(94.3, 15.7, 62.9, 126, 5.24,
                       31.4, 1.05, 1.05, 2.1, 10.5)) 
pumpData <- list(x = c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22)) 
pumpInits <- list(alpha = 1, beta = 1,
                  theta = rep(0.1, pumpConsts$N)) 
newPump <- nimbleModel(code = pumpCode, name = 'pump', constants = pumpConsts,
                    data = pumpData, inits = pumpInits) 

pumpMCEM <- buildMCEM(model = newPump,
                      latentNodes = 'theta',
                      burnIn = 100,
                      mcmcControl = list(adaptInterval = 20),
                      boxConstraints = list( list( c('alpha', 'beta'), 
                          limits = c(0, Inf) ) ), 
                      buffer = 1e-6)

pumpMCEM(maxit = 20, m1 = 250, m2 = 500)
pumpMCEM(maxit = 50, m1 = 1000, m2 = 5000)


## testing adding numeric(), integer(), array(), matrix()
library(nimble)
Rfun <- nimbleFunction(run = function() {
    ##ans <- numeric(10, value = 2)
    ##
    #x <- 100 
    #ans <- integer(x) 
    #for(i in 1:x) {
    #    ans[i] <- i 
    #}
    ##
    ##ans <- matrix(1, nrow = 10, ncol = 1)
    ##y <- 3
    ##ans <- numeric(10, value = y) 
    ##ans <- array(y, dim = 10) 
    ##ans <- array(y, dim = c(10))
    ##z <- numeric(10)
    ##z[5] <- 3
    ##x <- 20
    ##y <- 30
    ##ans <- integer(z[5], value = x + y) 
    ##ans <- array(x+y, dim = z[5], type = 'integer')
    ##ans <- array(x+y, dim = c(z[5]), type = 'integer')
    ##
    ##x <- array(0, c(4,5))
    ##ans <- matrix(0, nrow = dim(x)[1], ncol = dim(x)[2]) 
    ##ans <- array(0, dim = c(dim(x)[1], dim(x)[2]))
    ##
    x <- 1
    y <- 2
    z <- 3
    ans <- array(0, dim = c(x, y, z))
    ##returnType(double(1))
    ##returnType(integer(1))
    ##returnType(double(2))
    returnType(double(3))
    return(ans)
})
Cfun <- compileNimble(Rfun)

Rfun()
Cfun()
class(Rfun()[1])
class(Cfun()[1])
class(Rfun()[1,1])
class(Cfun()[1,1])
class(Rfun()[1,1,1])
class(Cfun()[1,1,1])

## creates a length-100 integer vector, containing 1, 2, ..., 100 

## creates a 10x1 ones-matrix 

## the following three lines are equivalent 
## each creates a length-10 vector, with all elements equal to y 

## the following two lines are equivalent 
## each creates an integer vector of length z[5], with all elements equal to x+y 

## the following two lines are equivalent 
## each one creates a matrix of 0's of the same size as matrix x 

## the following creates a 3-dimensional array of 0's 



## testing adding numeric(), integer(), array(), matrix()
library(nimble)
library(testthat)
##source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

expected <- numeric(10)
Rfun <- nimbleFunction(run = function() {
    ans <- numeric(10)
    returnType(double(1))
    return(ans)
})
Cfun <- compileNimble(Rfun)
test_that('numeric', expect_equal(Rfun(), expected))
test_that('numeric', expect_equal(Cfun(), expected))
test_that('numeric', expect_identical(class(Rfun()[1]), 'numeric'))
test_that('numeric', expect_identical(class(Cfun()[1]), 'numeric'))


expected <- rep(3, length = 2)
Rfun <- nimbleFunction(run = function() {
    ans <- numeric(value = 3, length = 2)
    returnType(double(1))
    return(ans)
})
Cfun <- compileNimble(Rfun)
test_that('numeric', expect_equal(Rfun(), expected))
test_that('numeric', expect_equal(Cfun(), expected))
test_that('numeric', expect_identical(class(Rfun()[1]), 'numeric'))
test_that('numeric', expect_identical(class(Cfun()[1]), 'numeric'))

expected <- rep(9, 3)
Rfun <- nimbleFunction(run = function() {
    x <- numeric(10, value = 3)
    ans <- integer(x[2], x[2]+x[3]*2)
    returnType(integer(1))
    return(ans)
})
Cfun <- compileNimble(Rfun)
test_that('integer', expect_equal(Rfun(), expected))
test_that('integer', expect_equal(Cfun(), expected))
test_that('integer', expect_identical(class(Rfun()[1]), 'integer'))
test_that('integer', expect_identical(class(Cfun()[1]), 'integer'))

expected <- array(4, c(10,11))
Rfun <- nimbleFunction(run = function() {
    ans <- array(4, c(10,11))
    returnType(double(2))
    return(ans)
})
Cfun <- compileNimble(Rfun)
test_that('integer', expect_equal(Rfun(), expected))
test_that('integer', expect_equal(Cfun(), expected))
test_that('integer', expect_identical(class(Rfun()[1]), 'numeric'))
test_that('integer', expect_identical(class(Cfun()[1]), 'numeric'))

expected <- matrix(as.integer(0), nrow=4, ncol=5)
Rfun <- nimbleFunction(run = function() {
    x <- 4
    y <- 5
    ans <- matrix(init=FALSE, nrow=x, ncol=y, type='integer')
    returnType(integer(2))
    return(ans)
})
Cfun <- compileNimble(Rfun)
test_that('integer', expect_equal(Rfun(), expected))
test_that('integer', expect_equal(Cfun(), expected))
test_that('integer', expect_identical(class(Rfun()[1,1]), 'integer'))
test_that('integer', expect_identical(class(Cfun()[1,1]), 'integer'))




expected
Rfun()
Cfun()


## testing adding numeric(), integer(), array(), matrix()
library(nimble)
nfDef <- nimbleFunction(
    setup = function() {},
    run = function() {
        x <- 2
        bbb <- array(3, x+10, type='integer')
        ##print(bbb)
        bbb2 <- array(value=3, dim=c(3,4), type='double')
        ##print(bbb2)
        ccc <- array(value=3, dim=c(3,4), type='double', init=FALSE)
        ##print(ccc)
        x <- 5
        y <- matrix(2, ncol = 3, nrow = x)
        ##print(y)
        zxc <- array(2, c(x, dim(y)[2], y[2,2]))
        zxc[1,1,1] <- 98
        ##returnType(double(3));  return(zxc)
        zz <- numeric(4)
        zz[1] <- 99
        ##print(zz)
        z <- numeric(value=.5, 4)
        ##print(z)
        zz[2:4] <- z[1:3] + numeric(3, 10)
        ##print(zz)
        mat <- matrix(0, 5, 5)
        mat[1, 1] <- 97
        ##print(mat)
        arr <- array(4, type = 'integer', init = FALSE, dim=c(3,4))
        ##print(arr)
        ar2 <- array(dim = dim(z)[1])
        ##print(ar2)
        ar3 <- array(dim = c(dim(z)[1]))
        ##print(ar3)
        ar4 <- array(dim = c(dim(arr)[1], dim(arr)[2]), type = 'integer')
        ##print(ar4)
    }
)
Rnf <- nfDef()
##Rnf$run
Cnf <- compileNimble(Rnf, dirName)

Rnf$run()

Cnf$run()



tempdir()

##Cnf$run()



## ???

library(nimble)
sampler_dt <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        browser()
        1
        2
        3
        ###  node list generation  ###
        calcNodes  <- model$getDependencies(target)
    },
    run = function() {
        simulate(model, target)
        calculate(model, calcNodes)
        nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    },
    methods = list(
        reset = function() { }
    ), where = getLoadingNamespace()
)
pumpCode <- nimbleCode({
	for(i in 1:N){
		theta[i] ~ dgamma(alpha,beta)
		lambda[i] <- theta[i]*t[i]
		x[i] ~ dpois(lambda[i])
	}
	alpha ~ dexp(1.0)
        alpha2 ~ dnorm(alpha, 1)
	beta ~ dpois(lambda=3)  # Make beta discrete valued
})
pumpConsts <- list(N=10,t=c(94.3,15.7,62.9,126,5.24,31.4,1.05,1.05,2.1,10.5))
pumpData <- list(x=c(5,1,5,14,3,19,1,1,4,22))
pumpInits <- list(alpha=1,beta=1,theta=rep(0.1,pumpConsts$N))
Rmodel <- nimbleModel(code=pumpCode,name='pump',constants=pumpConsts,data=pumpData,inits=pumpInits)

spec <- configureMCMC(Rmodel)
spec$getMonitors()
spec$addMonitors(c('theta', 'alpha2'))
spec$getMonitors()
spec$printSamplers()
spec$addSampler(type = 'dt', target = 'alpha')
spec$printSamplers()

Rmcmc <- buildMCMC(spec)

## testing new function in modelBaseClass: model$getNodeFunctions()
library(nimble)
pumpCode <- nimbleCode({
	for(i in 1:N){
		theta[i] ~ dgamma(alpha,beta)
		lambda[i] <- theta[i]*t[i]
		x[i] ~ dpois(lambda[i])
	}
	alpha ~ dexp(1.0)
        alpha2 ~ dnorm(alpha, 1)
	beta ~ dpois(lambda=3)  # Make beta discrete valued
})
pumpConsts <- list(N=10,t=c(94.3,15.7,62.9,126,5.24,31.4,1.05,1.05,2.1,10.5))
pumpData <- list(x=c(5,1,5,14,3,19,1,1,4,22))
pumpInits <- list(alpha=1,beta=1,theta=rep(0.1,pumpConsts$N))
Rmodel <- nimbleModel(code=pumpCode,name='pump',constants=pumpConsts,data=pumpData,inits=pumpInits)
Cmodel <- compileNimble(Rmodel)

nodes <- c('theta[3]')
nodes <- c('theta[3]', 'x[10]')
nodes <- c('theta[3]', 'x[9:10]')
Rmodel$getNodeFunctions(nodes)$calculate
Cmodel$getNodeFunctions(nodes)[[1]]$calculate

Rmodel$getNodeNames()
modelDef <- Rmodel$getModelDef()
gids <- mdef$nodeName2GraphIDs(nodes)
gids
mdef$graphIDs2indexedNodeInfo(gids)
Rmodel$nodeFunctions[[2]]$calculate


## making dmnorm conjugate sampler work with dependent
## nodes of different size from target node
library(nimble)
source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

code <- nimbleCode({
    x[1:3] ~ dmnorm(mu0[1:3], prec = ident[1:3,1:3])
    mu_y2[1:2] <- asCol(a[1:2]) + B[1:2,1:3] %*% asCol(x[1:3])
    mu_y3[1:3] <- asCol(a[1:3]) + B[1:3,1:3] %*% asCol(x[1:3])
    mu_y5[1:5] <- asCol(a[1:5]) + B[1:5,1:3] %*% asCol(x[1:3])
    y2[1:2] ~ dmnorm(mu_y2[1:2], prec = prec_y[1:2,1:2])
    y3[1:3] ~ dmnorm(mu_y3[1:3], prec = prec_y[1:3,1:3])
    y5[1:5] ~ dmnorm(mu_y5[1:5], prec = prec_y[1:5,1:5])
})

mu0 <- rep(0,3)
ident <- diag(3)
a <- 11:15
B <- matrix(1:15, nrow=5, ncol=3, byrow=TRUE)
prec_y <- diag(1:5)

constants <- list(mu0=mu0, ident=ident, a=a, B=B, prec_y=prec_y)
data <- list(y2=1:2, y3=1:3, y5=1:5)
inits <- list(x=rep(0,3))

Rmodel <- nimbleModel(code, constants, data, inits)

spec <- configureMCMC(Rmodel)
##spec$getSamplers()
Rmcmc <- buildMCMC(spec)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Rmcmc$run(10)

set.seed(0)
Cmcmc$run(10)

Rsamples <- as.matrix(Rmcmc$mvSamples)
Csamples <- as.matrix(Cmcmc$mvSamples)

Rsamples - Csamples

Rsamples[c(1,10),]
Csamples[c(1,10),]
##         x[1]       x[2]      x[3]
##[1,] 4.966869 -1.1252473 -3.922769
##[2,] 4.972376 -0.5288582 -4.511282


##node <- quote(mu0[1:3])
##node <- quote(prec_y[1:2, 1:2])
##node <- quote(prec_y[1:3, 1:3])
##node <- quote(prec_y[1:5, 1:5])
##nimble:::cc_expandDetermNodesInExpr(Rmodel, node)
##debug(nimble:::cc_expandDetermNodesInExpr)
##undebug(nimble:::cc_expandDetermNodesInExpr)
##debug(Rmcmc$samplerFunctions$contentsList[[1]]$run)
##a <- spec$getSamplerDefinition(1)
##createNamedObjectsFromList(a, writeToFile = '~/temp/del.R')







## ???
library(nimble)

d <- read.csv("http://personal.bgsu.edu/~albert/data/gateway.csv")
library(dplyr)
d$month <- as.numeric(d$Month)
d$year <- d$Year - 2002

worshipCode <- nimbleCode({ 
  for (i in 1:N){
  y[i] ~ dpois(lambda[i])
  log(lambda[i]) <- mu + epsilon[i] + alpha[month[i]]
  epsilon[i] ~ dnorm(0, sd=sigma.y)
  }
  for (j in 1:J){
  alpha[j] ~ dnorm(0, sd=sigma.month)
}
mu ~ dnorm(0, sd=1000)
sigma.y ~ dunif(0, 100)
sigma.month ~ dunif(0, 100)
})

worshipConsts <- list(N=484,
                      J=12,
                      K=10,
                      month=d$month)
worshipData <- list(y=d$Count)
worshipInits <- list(mu=5, sigma.y=.4, sigma.month=.4,
                     alpha=rep(0, 12))

worship <- nimbleModel(code = worshipCode, name = 'worship', 
                    constants = worshipConsts,
                    data = worshipData,
                    inits=worshipInits)

Rmodel <- worship

spec <- configureMCMC(Rmodel)
spec$printSamplers()
Rmcmc <- buildMCMC(spec)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Cmcmc$run(10000)
samples <- as.matrix(Cmcmc$mvSamples)
apply(samples, 2, mean)


## error produced in model$check() when using equals()
library(nimble)
code <- nimbleCode({
    for (i in 1:5) {
        x[i] ~ dbern(0.5)
        ##y[i] <- equals(x[i], 0)
    }
    sss <- sum(x[1:5])
    none <- equals(sss, 0)
})
constants <- list()
data <- list()
inits <- list(x = c(0,0,1,1,1))

Rmodel <- nimbleModel(code, constants, data, inits)

Rmodel$x
Rmodel$sss
Rmodel$none

## testing new log=TRUE option of RW sampler

library(nimble)
code <- nimbleCode({
    sigma ~ dunif(0,10)
    a ~ dnorm(0, sd=sigma)
    b ~ dnorm(sigma, tau=3)
})
constants <- list()
data <- list(a=1, b=4)
inits <- list(sigma=1)

Rmodel <- nimbleModel(code, constants, data, inits)
Cmodel <- compileNimble(Rmodel)
spec <- configureMCMC(Rmodel, nodes=NULL)

spec$addSampler('sigma', 'RW_log')
### --- or ---
spec$addSampler('sigma', 'RW', control=list(log=TRUE))

Rmcmc <- buildMCMC(spec)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Rmcmc$run(10)
as.matrix(Rmcmc$mvSamples)

set.seed(0)
Cmcmc$run(10)
as.matrix(Cmcmc$mvSamples)

##         sigma
## [1,] 3.535852
## [2,] 3.535852
## [3,] 3.535852
## [4,] 5.352671
## [5,] 5.352671
## [6,] 3.986347
## [7,] 3.963423
## [8,] 3.963423
## [9,] 3.963423
##[10,] 3.963423




## Github issue #107
## initializeModel issues false warning about RHSonly variable not initialized
## reduced case
library(nimble)

code <- nimbleCode({
    for (n in 1:N){
        a0[n,1:3] ~ dmnorm(mu[1:3], Sigma[1:3,1:3])
        ##lambda[n,1:3] <- exp(a0[n,1:3]) ## alternative, which shifts the problem to lambda
        for (k in 1:K){
            lambda[n,k] <- exp(a0[n,k])
            y[n,k]~dpois(lambda[n,k])
        }
    }
})

constants <- list(N = 3, K = 3)
data <- list(y = matrix(rpois(9, 3), nrow=3))
inits <- list(mu = rep(0,3), Sigma = diag(3))

Rmodel <- nimbleModel(code, constants, data, inits)

Cmodel <- compileNimble(Rmodel)

Rmcmc <- buildMCMC(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

Rmcmc$run(1)
Cmcmc$run(1)



## testing new binary sampler
library(nimble)

code <- nimbleCode({
    a ~ dbern(0.5)
    b ~ dbern(0.6)
    c ~ dbern(0.05)
    d ~ dbin(prob=0.2, size=1)
    e ~ dbinom(prob=0.9, size=1)
    f ~ dbern(0.5)
    g ~ dbern(0.5)
    h ~ dbern(0.5)
    for(i in 1:10)
        yf[i] ~ dnorm(f, sd = 1)
    for(i in 1:10)
        yg[i] ~ dnorm(g, sd = 1)
    for(i in 1:10)
        yh[i] ~ dnorm(h, sd = 1)
    ##x ~ dnorm(0,1)
    ##y ~ dnorm(x*x, 1)
    ##z ~ dnorm(x*y, 1)
    ##zz ~ dnorm(z + z, 1)
})
constants <- list()
data <- list(yf = c(rep(0,2), rep(1,8)), yg = c(rep(0,8), rep(1,2)), yh = c(rep(0,5), rep(1,5)))
inits <- list(a=0, b=0, c=0, d=0, e=0, f=0, g=0, h=0)
Rmodel <- nimbleModel(code, constants, data, inits)

##spec <- configureMCMC(Rmodel, autoBlock = TRUE)
##spec$printSamplers()

Rmodel$isBinary('a')
Rmodel$isBinary('b')
Rmodel$isBinary('c')
Rmodel$isBinary('d')
Rmodel$isBinary('e')
Rmodel$isBinary('f')
Rmodel$isBinary('g')
Rmodel$isBinary('h')

spec <- configureMCMC(Rmodel, nodes = NULL)
spec$addSampler('a', 'binary', print=FALSE)
spec$addSampler('b', 'binary', print=FALSE)
spec$addSampler('c', 'binary', print=FALSE)
spec$addSampler('d', 'binary', print=FALSE)
spec$addSampler('e', 'binary', print=FALSE)
spec$addSampler('f', 'binary', print=FALSE)
spec$addSampler('g', 'binary', print=FALSE)
spec$addSampler('h', 'binary', print=FALSE)
spec$printSamplers()
Rmcmc <- buildMCMC(spec)
##Rmcmc$run(3)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Cmcmc$run(100000)
samples <- as.matrix(Cmcmc$mvSamples)
means <- apply(samples, 2, mean)
means


# Slice sampler name conflict with library(bbmle) 'slice' function
library(bbmle)
library(nimble)  # No messages about slice are generated

pumpCode <- nimbleCode({
	for(i in 1:N){
		theta[i] ~ dgamma(alpha,beta)
		lambda[i] <- theta[i]*t[i]
		x[i] ~ dpois(lambda[i])
	}
	alpha ~ dexp(1.0)
        alpha2 ~ dnorm(alpha, 1)
	beta ~ dpois(lambda=3)  # Make beta discrete valued
})

pumpConsts <- list(N=10,t=c(94.3,15.7,62.9,126,5.24,31.4,1.05,1.05,2.1,10.5))
pumpData <- list(x=c(5,1,5,14,3,19,1,1,4,22))
pumpInits <- list(alpha=1,beta=1,theta=rep(0.1,pumpConsts$N))
Rmodel <- nimbleModel(code=pumpCode,name='pump',constants=pumpConsts,data=pumpData,inits=pumpInits)

spec <- configureMCMC(Rmodel)

spec$getSamplers()
Rmcmc <- buildMCMC(spec)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Cmcmc$run(2000000)
samples <- as.matrix(Cmcmc$mvSamples)
apply(samples, 2, mean)


## testing passing a compiled model as a setup argument to a NF

library(nimble)

code <- nimbleCode({
     a ~ dnorm(0, 1)
})
constants <- list()
data <- list()
inits <- list(a = 0)

Rmodel <- nimbleModel(code, constants, data, inits)
Cmodel <- compileNimble(Rmodel)

nfDef <- nimbleFunction(
    setup = function(model, node) {},
    run = function() {
        simulate(model, node)
    }
)

Rnf <- nfDef(Rmodel, 'a')
Rnf <- nfDef(Cmodel, 'a')

Cnf <- compileNimble(Rnf, project = Rmodel)

set.seed(0)
Rmodel$a
Rnf$run()
Rmodel$a

set.seed(0)
Cmodel$a
Cnf$run()
Cmodel$a



## testing passing a compiled model object to configureMCMC(), or buildMCMC()

library(nimble)

code <- nimbleCode({
     a ~ dnorm(0, 1)
})
constants <- list()
data <- list()
inits <- list(a = 0)

Rmodel <- nimbleModel(code, constants, data, inits)
Cmodel <- compileNimble(Rmodel)

###
Rmcmc <- buildMCMC(Cmodel)

###
spec <- configureMCMC(Cmodel)
spec$getSamplers()
Rmcmc <- buildMCMC(spec)

Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Rmcmc$run(10)
as.matrix(Rmcmc$mvSamples)

set.seed(0)
Cmcmc$run(10)
as.matrix(Cmcmc$mvSamples)


## testing problem with covariance matrix in BUGS model causing a chol() error at time of model checking

library(nimble)

code <- nimbleCode({
    C[1,1] <- 1
    C[1,2] <- 0
    C[2,1] <- 0
    C[2,2] <- 1
    mu[1] <- 0
    mu[2] <- 0
    y[1:2] ~ dmnorm(mu[1:2], prec = C[1:2, 1:2])
})
constants <- list()
data <- list(y = c(0,0))
inits <- list()

Rmodel <- nimbleModel(code, constants, data, inits)





##If you're interested, here's the ecological/scientific background to the project: the potato psyllid spreads a disease called Zebra Chip Disease that is damaging potato and tomato crops in California and western US and Mexico. The disease is new, having only been described about 20 years ago, but the potato psyllid is native to the area. So what's changed to cause the disease outbreaks? I'm using data from where and when potato psyllid museum specimens were collected in California to test if psyllid populations have increased over the last century and if climate change may play a role. If so, this could explain the disease outbreaks. It's complicated though because museum data is messy, specifically, the museum specimens give me presence-only information. So I'm using lists of related species collected at the same time as the psyllids to infer psyllid absences, based on an assumed correlation between number of species collected and collecting effort. This is the "list_length" covariate that's in the model.

##BLOCK SAMPLE PARAMS THAT GO IN THE SAME LINEAR PREDICTOR,
##i.e. all the betapN terms in this one

{
    ## Priors
    ## For the site random effect
    for(j in 1:nsite) { 
        alpha[j] ~ dnorm(mu.alpha, tau.alpha) 
    }
    mu.alpha ~ dnorm(0, 0.001)
    tau.alpha <- 1 / (sigma.alpha * sigma.alpha)
    sigma.alpha ~ dunif(0, 5)
    
    ## Grand mean
    muq ~ dnorm(0, 0.001)
    
    ## For the fixed effect coefficients
    betaq ~ dnorm(0, 0.001)
    betap1 ~ dnorm(0, 0.001)
    betap2 ~ dnorm(0, 0.001)
    betap3 ~ dnorm(0, 0.001)
    betap4 ~ dnorm(0, 0.001)
    betap5 ~ dnorm(0, 0.001)
    betap6 ~ dnorm(0, 0.001)
    betap7 ~ dnorm(0, 0.001)
    betap8 ~ dnorm(0, 0.001)

    ## Likelihood
    for (i in 1:nlist){ # i = events (year-months)
        for(j in 1:nsite) { # j = sites

            logit(p[i,j]) <- betap1*year[i,j] + betap2*pow(year[i,j],2) + # Year quadratic effects
                betap2*month[i,j] + betap3*pow(month[i,j],2) + # month quadratic effects
                    betap4*aet[i,j] + betap5*cwd[i,j] + betap6*tmn[i,j] + betap7*tmx[i,j] + # Climate effects
                        alpha[j] # Random effects

            Y[i,j] ~ dbern(p[i,j]) # Occupancy probability

            logit(q[i,j]) <- Y[i,j] + muq + betaq*list_length[i,j] # Logistic regression for detection

            detectionMatrix[i,j] ~ dbern(q[i,j])  # Distribution for random part; observed presences relating to detection probability

        } #j
    } #i
}



## bug with forwardsolve in nimbleFunctions

library(nimble)
A <- matrix(c(2,1,0,3), nrow = 2)
b <- c(1,5)

nf <- nimbleFunction(
    run = function(A = double(2), b = double(1)) {
        declare(x, double(1))
        setSize(x, 2)
        x[1:2] <- forwardsolve(A[1:2, 1:2], b[1:2])
        ##x <- forwardsolve(A[1:2,1:2], b[1:2])
        x[1:2,1:2] <- chol(A[1:2,1:2])
        returnType(double(1))
        ##declare(x, double(2))
        ##setSize(x, 2, 2)
        ##returnType(double(2))
        return(x)
    })

cnf <- compileNimble(nf, dirName = '.')

nf(A,b)
cnf(A,b)



## bug with forwardsolve in BUGS models

library(nimble)

code <- nimbleCode({
    A[1,1] <- 2
    A[2,1] <- 1
    A[1,2] <- 0
    A[2,2] <- 3
    b[1] <- 1
    b[2] <- 5
    x[1:2] <- forwardsolve(A[1:2, 1:2], b[1:2])
})

Rmodel <- nimbleModel(code)
Cmodel <- compileNimble(Rmodel)

Cmodel$A
Cmodel$b
Cmodel$x
calculate(Cmodel)
## Process R floating point exception: 8 at Thu Apr  7 14:38:02 2016



## testing using forwardsolve and backsolve in a user-defined function
## in a BUGS model

library(nimble)
A <- array(NA, c(2,2))
A[1,1] <- 2
A[2,1] <- 1
A[1,2] <- 0
A[2,2] <- 3
b <- c(1,5)

myFS <- nimbleFunction(
    run = function(A = double(2), b = double(1)) {
        returnType(double(1))
        ret <- forwardsolve(A, b)
        return(ret)
    }
)

myChol <- nimbleFunction(
    run = function(A = double(2)) {
        returnType(double(2))
        ret <- chol(A)
        return(ret)
    }
)

myInverse <- nimbleFunction(
    run = function(A = double(2)) {
        returnType(double(2))
        ret <- inverse(A)
        return(ret)
    }
)

code <- nimbleCode({
    x1[1:2] <- forwardsolve(A[1:2, 1:2], b[1:2])  ## forwardsolve() directly here does NOT work
    x2[1:2] <- myFS(A = A[1:2, 1:2], b = b[1:2])
    C1[1:2,1:2] <- chol(A[1:2, 1:2])
    C2[1:2,1:2] <- myChol(A[1:2, 1:2])
    I1[1:2,1:2] <- inverse(A[1:2, 1:2])
    I2[1:2,1:2] <- myInverse(A[1:2, 1:2])
    xx1[1:2,1:2] <- chol(A[1:2, 1:2]) + 2 * inverse(A[1:2, 1:2])
    xx2[1:2,1:2] <- myInverse(A[1:2, 1:2]) * 2 + myChol(A[1:2, 1:2])
})

Rmodel <- nimbleModel(code, constants = list(A=A, b=b))
Cmodel <- compileNimble(Rmodel)

m <- Rmodel
m <- Cmodel

calculate(m)
m$x1; m$x2
m$C1; m$C2
m$I1; m$I2
m$xx1; m$xx2

## another example of MVN conjugate sampler, for test-mcmc.R
## using both cov and prec parametrizaions of MVN,
## and various linear links

library(nimble)

set.seed(0)
mu0 <- rep(0,5)
ident <- diag(5)
a <- array(rnorm(20), c(4,5))
B <- array(NA, c(4, 5, 5))
for(i in c(1,3)) {
    A <- array(rnorm(25,i), c(5,5))
    A <- A + t(A) + 5*i*diag(5)
    B[i,,] <- A
}
B[2,,] <- diag(5)
B[4,,] <- diag(5)
M_y <- array(NA, c(4, 5, 5))
for(i in 1:4)
    M_y[i,,] <- i * diag(5)
x <- array(0, c(4,5))
y <- array(0, c(4,5))

code <- nimbleCode({
    for(i in 1:2)
        x[i,1:5] ~ dmnorm(mu0[1:5], prec = ident[1:5,1:5])
    for(i in 3:4)
        x[i,1:5] ~ dmnorm(mu0[1:5], cov  = ident[1:5,1:5])
    for(i in 1:4)
        mu_y[i,1:5] <- asCol(a[i,1:5]) + B[i,1:5,1:5] %*% asCol(x[i,1:5])
    y[1,1:5] ~ dmnorm(mu_y[1,1:5], prec = M_y[1,1:5,1:5])
    y[2,1:5] ~ dmnorm(mu_y[2,1:5], cov  = M_y[2,1:5,1:5])
    y[3,1:5] ~ dmnorm(mu_y[3,1:5], prec = M_y[3,1:5,1:5])
    y[4,1:5] ~ dmnorm(mu_y[4,1:5], cov  = M_y[4,1:5,1:5])
})
constants <- list(mu0=mu0, ident=ident, a=a, B=B, M_y=M_y)
data <- list(y=y)
inits <- list(x=x)
Rmodel <- nimbleModel(code, constants, data, inits)
spec <- configureMCMC(Rmodel)
spec$getSamplers()
Rmcmc <- buildMCMC(spec)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Rmcmc$run(10)
Rsamples <- as.matrix(Rmcmc$mvSamples)

set.seed(0)
Cmcmc$run(10)
Csamples <- as.matrix(Cmcmc$mvSamples)

Rsamples - Csamples

Rsamples[c(1,10),]
Csamples[c(1,10),]
##        x[1, 1]     x[2, 1]     x[3, 1]    x[4, 1]     x[1, 2]   x[2, 2]
##[1,] -0.7098220 -0.92579244  0.02859492 -0.9271968 -0.03440439 0.1133790
##[2,] -0.6642612  0.09915454 -0.13079585  0.6298495 -0.08432850 0.8806305
##        x[3, 2]    x[4, 2]     x[1, 3]   x[2, 3]     x[3, 3]    x[4, 3]
##[1,] 0.09213596  0.3257091  0.27465408 -1.410726 -0.12937422 -0.2311789
##[2,] 0.01777368 -0.6051182 -0.09322287 -1.723685  0.07759651  0.4119755
##       x[1, 4]     x[2, 4]     x[3, 4]     x[4, 4]   x[1, 5]  x[2, 5]
##[1,] 0.5071505 0.140652780  0.01778143  0.51383003 0.3639871 2.146307
##[2,] 0.4781541 0.002361055 -0.02517465 -0.06995162 0.4202336 1.225317
##         x[3, 5]    x[4, 5]
##[1,] -0.02913630 -0.2563026
##[2,]  0.04303507  0.2716980




## Richard McElreath Statistical Rethinking package,
## might have lots of good data and hierarchical (Bayesian) modeling examples
install.packages('rethinking')
library(rethinking)



## tesing MVN conjugate sampler

library(nimble)

set.seed(0)
mu0 = 1:3
Q0 = matrix(c(1, .2, .8, .2, 2, 1, .8, 1, 2), nrow = 3)
Q = solve(matrix(c(3, 1.7, .9, 1.7, 2, .6, .9, .6, 1), nrow = 3))
a = c(-2, .5, 1)
B = matrix(rnorm(9), 3)

code <- nimbleCode({
  mu[1:3] ~ dmnorm(mu0[1:3], Q0[1:3, 1:3])
  y_mean[1:3] <- asCol(a[1:3]) + B[1:3, 1:3] %*% asCol(mu[1:3])
  y[1:3] ~ dmnorm(y_mean[1:3], Q[1:3, 1:3])
})
mu <- mu0 + chol(solve(Q0)) %*% rnorm(3)
y <- c(a + B%*%mu + chol(solve(Q)) %*% rnorm(3))
data = list(mu0 = mu0, Q0 = Q0, Q = Q, a = a, B = B, y = y)
muQtrue = t(B) %*% Q%*%B + Q0
muMeanTrue = c(solve(muQtrue, crossprod(B, Q%*%(y-a)) + Q0%*%mu0))

Rmodel <- nimbleModel(code, data= data)
spec <- configureMCMC(Rmodel)
spec$getSamplers()
##spec$getSamplerDefinition(1)
Rmcmc <- buildMCMC(spec)
Cmodel <- compileNimble(Rmodel)

Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Rmcmc$run(10)
Rsamples <- as.matrix(Rmcmc$mvSamples)

set.seed(0)
Cmcmc$run(10)
Csamples <- as.matrix(Cmcmc$mvSamples)

Rsamples
Csamples
##         mu[1]    mu[2]    mu[3]
## [1,] 3.877574 2.201065 1.084173
## [2,] 2.537805 1.916830 2.084021
## [3,] 4.426193 2.326982 1.567054
## [4,] 2.456652 1.921560 1.892777
## [5,] 2.959890 2.133226 1.506502
## [6,] 3.084360 1.560040 1.941621
## [7,] 3.178667 2.072212 2.611881
## [8,] 2.996492 2.210513 2.795391
## [9,] 2.494759 1.539610 2.118232
##[10,] 2.880285 1.826356 1.805385

Rsamples - Csamples
##              mu[1] mu[2] mu[3]
## [1,] -8.881784e-16     0     0
## [2,]  0.000000e+00     0     0
## [3,] -8.881784e-16     0     0
## [4,] -4.440892e-16     0     0
## [5,] -4.440892e-16     0     0
## [6,]  0.000000e+00     0     0
## [7,] -4.440892e-16     0     0
## [8,] -4.440892e-16     0     0
## [9,]  0.000000e+00     0     0
##[10,] -4.440892e-16     0     0


## inconsistent behavior between R and C chol()

n <- 4
set.seed(0)
tmp <- array(rnorm(n^2), c(n,n)) + n*diag(n)
A <- tmp + t(tmp)
R <- chol(A)
L <- t(R)
B <- array(as.numeric(1:n^2), c(n,n))
b <- as.numeric(1:n)
mat <- A + L

sym <- function(mat, lower) {
    ind <- if(lower) lower.tri(mat) else upper.tri(mat)
    ans <- mat
    ans[!ind] <- 0
    ans <- ans + t(ans) + diag(diag(mat))
    ans
}

## R execution
chol(mat)
chol(sym(mat, FALSE))
chol(mat) - chol(sym(mat, FALSE)) ## R takes the UPPER part of the matrix

## Rnf execution
Rnf$chol_oPA_plus_L_cP() - chol(sym(mat, FALSE))  ## Rnf also UPPER part of matrix

## Cnf execution
Cnf$chol_oPA_plus_L_cP() - chol(sym(mat, FALSE))  ## Cnf now also uses UPPER part of matrix





## test case for forwardsolve, backsolve in the DSL
## A.triangularView<Eigen::Upper>().solve(b);
## could make a file called nimble_macros.h:
##define ......

library(nimble)
## inputs
n <- 4
tests <- c(
    'chol(A)',
    'forwardsolve(A, b)', 'forwardsolve(A, B)', 'forwardsolve(L, b)', 'forwardsolve(L, B)',
       'backsolve(A, b)',    'backsolve(A, B)',    'backsolve(R, b)',    'backsolve(R, B)',
           'solve(A, b)',        'solve(A, B)',        'solve(R, b)',        'solve(R, B)'
)
## derived objects
testNames <- as.character(sapply(tests, nimble:::Rname2CppName))
testDims <- sapply(grepl('_b_', testNames), function(x) if(x) 1 else 2)
methods <- mapply(function(test, dim) eval(substitute(function() { ans <- TEST; returnType(double(DIM)); return(ans) }, list(TEST = parse(text=test)[[1]], DIM = dim))), tests, testDims)
names(methods) <- testNames
set.seed(0)
tmp <- array(rnorm(n^2), c(n,n)) + n*diag(n)
A <- tmp + t(tmp)
R <- chol(A)
L <- t(R)
B <- array(as.numeric(1:n^2), c(n,n))
b <- as.numeric(1:n)
nfDef <- nimbleFunction(
    setup = function(A, B, b) {
        R <- chol(A)
        L <- t(R)
    },
    run = function() {},
    methods = methods
)
Rnf <- nfDef(A, B, b)

Cnf <- compileNimble(Rnf)

testOneCase <- function(test, testName, Rnf, Cnf) {
    Rans <- eval(parse(text=test)[[1]])
    Rnfans <- eval(substitute(Rnf$TEST(), list(TEST=as.name(testName))))
    Cnfans <- eval(substitute(Cnf$TEST(), list(TEST=as.name(testName))))
    dif <- max(abs(Rans - Rnfans))
    if(dif > 1E-15) return(1)
    dif <- max(abs(Rans - Cnfans))
    if(dif > 1E-15) return(1)
    return(0)
}
for(i in seq_along(tests)) {
    test <- tests[i]
    testName <- testNames[i]
    testResult <- testOneCase(test, testName, Rnf, Cnf)
    if(testResult != 0) message('failed test: ', test)
}




library(nimble)

n <- 3
set.seed(0)
tmp <- array(rnorm(n^2), c(n,n)) + diag(n)
A <- tmp + t(tmp)
B <- array(1:n^2, c(n,n))
b <- as.numeric(1:n)


Rnf <- nimbleFunction(
    run = function(A = double(2), b = double(1), B = double(2)) {
        ##ans2 <- forwardsolve(A, b)
        ##ans3 <- forwardsolve(A, B)
        ans4 <- chol(A + A)
    }
)

Cnf <- compileNimble(Rnf)   ## this works




nfDef <- nimbleFunction(
    setup = function(A, B, b) {},
    run = function() {
        ans2 <- forwardsolve(A, b)
        ##ans3 <- forwardsolve(A, B)
    }
)

Rnf <- nfDef(A, B, b)
Cnf <- compileNimble(Rnf)   ## FAIL


library(nimble)

n <- 3
set.seed(0)
tmp <- array(rnorm(n^2), c(n,n)) + diag(n)
A <- tmp + t(tmp)
B <- array(1:n^2, c(n,n))
b <- as.numeric(1:n)


Rnf <- nimbleFunction(
    run = function(A = double(2), b = double(1), B = double(2)) {
        ans1 <- chol(A)
        print('chol(A)')
        print(ans1)
        ##
        ##ans2 <- forwardsolve(A, b)
        ##print('forwardsolve(A, b)')
        ##print(ans2)
        ####
        ##ans3 <- forwardsolve(A, B)
        ##print('forwardsolve(A, B)')
        ##print(ans3)
        ####
        ##ans4 <- backsolve(A, b)
        ##print('backsolve(A, b)')
        ##print(ans4)
        ####
        ans5 <- backsolve(A, B)
        print('backsolve(A, B)')
        print(ans5)
    }
)

Cnf <- compileNimble(Rnf)


set.seed(0)
tmp <- array(rnorm(9), c(3,3)) + diag(3)
A <- tmp + t(tmp)
B <- array(1:9, c(3,3))
b <- 1:3

chol(A)
forwardsolve(A, b)
forwardsolve(A, B)
backsolve(A, b)
backsolve(A, B)

Rnf(A, b, B)
Cnf(A, b, B)

##debug(exprClasses_labelForEigenization)
##undebug(exprClasses_labelForEigenization)
##debug(exprClasses_eigenize)
##undebug(exprClasses_eigenize)
##writeCode(nimDeparse(code))
##writeCode(nimDeparse(compileInfo$nimExpr))


## using compareMCMCs, and MCMC comparisons

library(nimble)

code1 <- nimbleCode({ ## really toy model
    a ~ dnorm(0, 1)
    b ~ dnorm(a, 1)
    sigma ~ dunif(0.5, 1)
    for(i in 1:3) y[i] ~ dnorm(b, sd = sigma)
})

input1 <- list(code = code1, data = list(y = 1:3), inits = list(a = 0.5)) ## can also provide constants
results1 <- compareMCMCs(input1, MCMCs = c('nimble')) ## run a single case

results2 <- compareMCMCs(input1, MCMCs = c('nimble','noConj')) ## run two more cases
results2[[1]]$efficiency ## inspect
results2[[1]]$summary
results2[[1]]$timing
make_MCMC_comparison_pages(results2, 'model1b') ## could generate a comparison of these

## can combine results as follows
## If the same name (“nimble” in this case) was used multiple times one will have to be renamed:
results2[[1]] <- rename_MCMC_comparison_method('nimble', 'another nimble', results2[[1]])
results3 <- combine_MCMC_comparison_results(results1[[1]], results2[[1]], name = 'combined results')
make_MCMC_comparison_pages(results3, 'model1c')



## testing calc_dmnormConjugacyContributions

nimble:::identityMatrix(3)

calc_dmnormConjugacyContributions(diag(3), diag(3), 1)
calc_dmnormConjugacyContributions(diag(3), 7*diag(3), 1)
calc_dmnormConjugacyContributions(diag(3), 8*diag(3), 2)
calc_dmnormConjugacyContributions(3*diag(3), 1+diag(3), 1)
calc_dmnormConjugacyContributions(3*diag(3), 2+diag(3), 2)

Cnf <- compileNimble(calc_dmnormConjugacyContributions)

Cnf(diag(3), diag(3), 1)
Cnf(diag(3), 7*diag(3), 1)
Cnf(diag(3), 8*diag(3), 2)
Cnf(3*diag(3), 1+diag(3), 1)
Cnf(3*diag(3), 2+diag(3), 2)



## compiler lesson 101
library(nimble)

nf <- nimbleFunction(
    setup = function() {
        a <- 1
    },
    run = function(arg = double()) {
        b <- a + arg
        returnType(double())
        return(b)
    })

Rnf <- nf()

## nfProc object is created every time we compile a new NF generator
## nfCompileInfo data structure for compilation info RCfunctionCompileClass"
## RCfunProcessing class to manage compilation info for a *single* run/member function
proj <- nimble:::nimbleProjectClass()
nfProc <- nimble:::nfProcessing(Rnf, project = proj)
nfProc$process(control = list(debug = TRUE))

## place to define 1st arg switching at top of genCpp_sizeProcessing
## switch name to nimArr_xxx in genCpp_processSpecificCalls






## testing of new parsed posterior text in conjugacy system
posteriorText = '{a<-1;
b
c
dnorm(mean = (prior_mean*prior_tau + contribution_mean) / (prior_tau + contribution_tau),
                            sd   = (prior_tau + contribution_tau)^(-0.5))
}'

posteriorText = '{ R <- chol(prior_prec + contribution_prec)
                        A <- prior_prec %*% asCol(prior_mean) + asCol(contribution_mean)
                        mu <- backsolve(R, forwardsolve(t(R), A))[,1]
                        dmnorm_chol(mean = mu, cholesky = R, prec_param = 1) }'

parsedTotalPosterior <- parse(text = posteriorText)[[1]]
parsedTotalPosterior
if(parsedTotalPosterior[[1]] != '{') parsedTotalPosterior <- substitute({POST}, list(POST = parsedTotalPosterior))
parsedTotalPosterior
prePosteriorCodeBlock <- parsedTotalPosterior[-length(parsedTotalPosterior)]
prePosteriorCodeBlock
posteriorExpr <- parsedTotalPosterior[[length(parsedTotalPosterior)]]
posteriorExpr
rDistribution <- cc_makeRDistributionName(as.character(posteriorExpr[[1]]))
rDistribution
dDistribution <- as.character(posteriorExpr[[1]])
dDistribution
argumentExprs <- as.list(posteriorExpr)[-1]
argumentExprs
argumentNames <- names(argumentExprs)
argumentNames
rCallExpr <- as.call(c(as.name(rDistribution), 1, argumentExprs))
rCallExpr
dCallExpr <- as.call(c(as.name(dDistribution), quote(VALUE), argumentExprs, log = 1))
dCallExpr
posteriorVars <- all.vars(parsedTotalPosterior)
posteriorVars
neededPriorParams <- gsub('^prior_', '', posteriorVars[grepl('^prior_', posteriorVars)])
neededPriorParams
neededContributionNames <- posteriorVars[grepl('^contribution_', posteriorVars)]
neededContributionNames
neededContributionDims <- inferContributionTermDimensions(prior)


## cases of dmnorm() parametrizations, and possible getParam() calls
library(nimble)

code <- nimbleCode({ x ~ dnorm(mu, var = v) })
inits <- list(x=0, mu=0, v=1)

set.seed(0)
tmp <- array(rnorm(9), c(3,3)) + diag(3)
A <- tmp + t(tmp)
code <- nimbleCode({ x[1:3] ~ dmnorm(mu[1:3], prec = Q[1:3, 1:3]) })
inits <- list(x=rep(0,3), mu=rep(0,3), Q=A)

md <- nimbleModel(code, inits = inits, returnDef = TRUE, debug=TRUE)

Rmodel <- nimbleModel(code, inits = inits)

Rmodel$getNodeNames()
nn <- Rmodel$getNodeNames(stochOnly = TRUE)
lifted <- Rmodel$getNodeNames(determOnly = TRUE)
nn
lifted
Rmodel$nodes[[lifted]]$simulate  ## calculates chol(Q) into lifted node
Rmodel$nodes[[nn]]$simulate
Rmodel$nodes[[nn]]$get_var
Rmodel$nodes[[nn]]$get_sd
Rmodel$nodes[[nn]]$get_tau
Rmodel$nodes[[nn]]$get_cholesky
Rmodel$nodes[[nn]]$get_prec
Rmodel$nodes[[nn]]$get_cov


## cases of dmnorm() parametrizations, and possible getParam() calls
library(nimble)

set.seed(0)
tmp <- array(rnorm(9), c(3,3)) + diag(3)
A <- tmp + t(tmp)

## parametrize in terms of prec matrix (Q):
code <- nimbleCode({
    x[1:3] ~ dmnorm(mu[1:3], prec = Q[1:3, 1:3])
})
inits <- list(x=rep(0,3), mu=rep(0,3), Q=A)

Rmodel <- nimbleModel(code, inits = inits)
Cmodel <- compileNimble(Rmodel)

Rmodel$getNodeNames()
nn <- Rmodel$getNodeNames(stochOnly = TRUE)
lifted <- Rmodel$getNodeNames(determOnly = TRUE)
nn
lifted
Rmodel[[lifted]]

Q <- A
ch <- chol(Q)
V <- solve(Q)
I <- diag(3)
Q
ch
V
t(ch) %*% ch
Q
backsolve(ch, forwardsolve(t(ch), I))
V

Rmodel$nodes[[nn]]$get_mean
Rmodel$nodes[[nn]]$get_mean()
Rmodel$getParam(nn, 'mean')

Rmodel$nodes[[nn]]$get_cholesky
Rmodel$nodes[[nn]]$get_cholesky()
Rmodel$getParam(nn, 'cholesky')
ch

Rmodel$nodes[[nn]]$get_prec
Rmodel$nodes[[nn]]$get_prec()
Rmodel$getParam(nn, 'prec')
Q

Rmodel$nodes[[nn]]$get_cov
Rmodel$nodes[[nn]]$get_cov()
Rmodel$getParam(nn, 'cov')
V

xx <- forwardsolve(t(ch), I)
t(xx) %*% xx

xx <- backsolve(ch, I)  ## would this work too???
xx %*% t(xx)

t(ch) %*% ch %*% V
t(ch) %*% ch
ch %*% V
Q


nfDef <- nimbleFunction(
    setup = function(model, node) {},
    run = function() {
        print('chol:')
        chh <- model$getParam(node, 'cholesky')
        print(chh)
        print('prec:')
        q <- model$getParam(node, 'prec')
        print(q)
        print('cov:')
        vv <- model$getParam(node, 'cov')
        print(vv)
    }
)

Rnf <- nfDef(Rmodel, nn)
Cnf <- compileNimble(Rnf, project=Rmodel)


ch
Q
V
Rnf$run()
Cnf$run()

## store something into model variable 'Q'
Rmodel$Q <- Q

## calculate deterministic dependents
Rmodel$nodes[[lifted]]$simulate  ## calculates chol(Q) into lifted node

Rmodel$nodes[[nn]]$simulate
Rmodel$nodes[[nn]]$get_cholesky
Rmodel$nodes[[nn]]$get_cholesky()
ch

## request 'prec'
Rmodel$nodes[[nn]]$get_prec  ## direct: model$Q
Rmodel$nodes[[nn]]$get_prec()
Q

## request 'cov'
## here's what we want to get cov: backsolve(ch,  forwardsolve(t(ch), I)  )
Rmodel$nodes[[nn]]$get_cov
Rmodel$nodes[[nn]]$get_cov()
V







    


## figuring out R / nimble functions:
## inverse, solve, forwardsolve, backsolve

library(nimble)
set.seed(0)
A <- array(rnorm(9), c(3,3))
A <- t(A) %*% A

nfDef <- nimbleFunction(
    setup <- function(A) {
        D <- diag(3)
        b <- 1:3
    },
    methods = list(
        inv = function() {
            Ainv <- inverse(A)
            print(Ainv)
        },
        ch = function() {
            R <- chol(A)
            print(R)
        },
        fs = function() {
            R <- chol(A)
            L <- t(R)
            Linv <- forwardsolve(L, D)
            print(Linv)
        },
        bs = function() {
            R <- chol(A)
        }
    )
)

Rnf <- nfDef(A)
Cnf <- compileNimble(Rnf)

print('inverse of A')
solve(A)
Rnf$inv()
Cnf$inv()

print('chol of A')
chol(A)
Rnf$ch()
Cnf$ch()

R <- chol(A)
L <- t(R)

print('forwardward solve')
forwardsolve(L, diag(3))

print('backward solve')
backsolve(R, diag(3))



## example of using nf$getDefinition() and nf_getDefinition()
library(nimble)
nfDef <- nimbleFunction(
    setup = function(a) { },
    run = function(b = double()) {
        a <<- a + b
        returnType(double())
        return(a)
    },
    methods = list(
        nothing = function() {
            print(a)
        })
)
Rnf <- nfDef(3)
Cnf <- compileNimble(Rnf)
getDefinition(nfDef)
Rnf$getDefinition()
Cnf$getDefinition()

## example of spec$getSamplerDefinition(..)
library(nimble)

code <- nimbleCode({
    a ~ dnorm(0, 1)
    b ~ dgamma(a, 1)
    for(i in 1:140)
        c[i] ~ dnorm(0,1)
})

Rmodel <- nimbleModel(code, check = FALSE)

spec <- configureMCMC(Rmodel, print = TRUE)
spec$printSamplers()
spec$getSamplerDefinition(1)
spec$getSamplerDefinition('b')

spec$getSamplers()
spec$printSamplers()
spec$addSampler(type='RW_block', target=c('a', 'b'))
spec$getSamplers()
spec$printSamplers()


## performance comparison of of new dynamic conjugate samplers, rats example
library(nimble)
rats <- readBUGSmodel('rats', dir = getBUGSexampleDir('rats'), returnModelComponentsOnly = TRUE)
Rmodel <- nimbleModel(rats$model, rats$data[c('N','T','x')], rats$data['Y'], rats$inits)
spec <- configureMCMC(Rmodel, print = TRUE)
out <- MCMCsuite(
    rats$model, rats$data[c('N','T','x')], rats$data['Y'], rats$inits,
    niter = 100000,
    MCMCs = c('staticConj', 'dynamicConj'),
    MCMCdefs = list(
        staticConj  = quote({ nimbleOptions(useDynamicConjugacy = FALSE)
                              configureMCMC(Rmodel) }),
        dynamicConj = quote({ nimbleOptions(useDynamicConjugacy = TRUE)
                              configureMCMC(Rmodel) })),
    makePlot = FALSE
)
sampdiff <- out$samples['staticConj',,] - out$samples['dynamicConj',,]
all(sampdiff == 0)
## [1] TRUE                      ## samples from each MCMC are identical
out$timing
## staticConj    dynamicConj
##     26.309         17.801     ## ~ 33% decrease in conjugate sampler runtime



## testing proper dmnorm conjugate calculations for Chris

## code = nimbleCode({
##    y[1:3]  ~ dmnorm(mu [1:3], cov = covy [1:3,1:3])
##    mu[1:3] ~ dmnorm(mu0[1:3], cov = covmu[1:3,1:3])
## })
## Concerned about some explicit inverses in conjugate update for mu,
## compared to use of already existing Cholesky factors.  In the R sampler,
## 1) we have: inverse(prior_prec + contribution_prec)
## 2) we use get_prec, which when we specify a dmnorm with a covariance
## has code like the following (in calc_dmnormAltParams()):
##  tmp <- t(cholesky) %*% cholesky
##  return(inverse(tmp))

library(nimble)

nimbleOptions('buildInterfacesForCompiledNestedNimbleFunctions')
nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)
nimbleOptions('buildInterfacesForCompiledNestedNimbleFunctions')
nimbleOptions()

code = nimbleCode({
   mu[1:3] ~ dmnorm(mu0[1:3], cov=covmu[1:3,1:3])
   y[1:3] ~ dmnorm(mu[1:3], cov = covy[1:3,1:3])
   a ~ dnorm(y[1], y[1])
   b ~ dgamma(a, y[1])
})

cov <- diag(3)
constants <- list(covy = cov, covmu = cov, mu0 = rep(0,3))
data <- list(y = 1:3)
inits <- list(mu = 4:6, a = 1, b = 1)

Rmodel <- nimbleModel(code, constants, data, inits)

spec <- configureMCMC(Rmodel)
spec$getSamplers()
spec$printSamplers()
spec$getSamplerDefinition('mu')
Rmcmc <- buildMCMC(spec)


Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

Rmcmc$samplerFunctions$contentsList[[2]]$scale
Cmcmc$samplerFunctions$contentsList[[2]]
Cmcmc$samplerFunctions$contentsList[[2]][['scale']]
Cmcmc$samplerFunctions$contentsList[[2]]$scale



set.seed(0)
Cmcmc$run(10000)
samples <- as.matrix(Cmcmc$mvSamples)
apply(samples, 2, mean)



## testing what works and what doesn't with nimbleFunctionLists, getParam, etc..

library(nimble)

getNimbleOption('checkModel')
nimbleOptions(checkModel = FALSE)
getNimbleOption('checkModel')

code <- quote({
    x ~ dgamma(1, 1)
    y[1] ~ dnorm(1, tau = 1*x)
    y[2] ~ dnorm(2, tau = 2*x)
    y[3] ~ dnorm(3, tau = 3*x)
})
Rmodel <- nimbleModel(code, data=list(), inits=list(x=1, y=11:13))

nfDef <- nimbleFunction(
    setup = function(model) {
        nl <- model$expandNodeNames('y')
        N <- length(nl)
        nfl <- nimbleFunctionList(node_stoch_dnorm)
        for (i in 1:N) nfl[[i]] <- model$nodeFunctions[[nl[i]]]
    },
    run = function() {
        for(i in 1:N) {
            print('i = ', i)
            val <- nfl[[i]]$get_mean()
            print('mean = ', val)
            val <- nfl[[i]]$get_tau()
            print('tau = ', val)
            val <- nfl[[i]]$get_var()
            print('var = ', val)
            val <- nfl[[i]]$get_value()
            print('value = ', val)
        }
        returnType(double(0))
        return(1)
    }
)
Rnf <- nfDef(Rmodel)
Rnf$run()
Cmodel <- compileNimble(Rmodel)
Cnf <- compileNimble(Rnf, project = Rmodel)
Cnf$run()


## testing of new dynamic conjugate samplers
### demo2 of check conjugacy
library(nimble)
code <- BUGScode({
    x ~ dbeta(3, 13)
    y[1] ~ dbin(x, 10)
    y[2] ~ dbin(x, 20)
    z ~ dbeta(10, 10)
    zz[1] ~ dbin(z, 10)
    zz[2] ~ dbin(z, 20)
    aaa ~ dbeta(10, 10)
    bbb[1] ~ dbin(aaa, 10)
    bbb[2] ~ dbin(aaa, 20)
    bbb[3] ~ dbin(aaa, 20)
    d ~ dnorm(0, 1)
    dd ~ dnorm(0, 1)
    ddd ~ dnorm(0, 1)
    e ~ dnorm(d + dd + ddd, 1)
    f ~ dlnorm(d + dd + ddd, 1)
    ff ~ dlnorm(d + dd, 1)
    g ~ dnorm(0, 1)
    gg ~ dnorm(5+g, 1)
})
constants <- list()
inits <- list(x = 0.5, z=0.6, aaa=0.3, d=0, dd=0, ddd=0, g=0)
data = list(y = c(3,4), zz=c(5,5), bbb=c(3,3,3), ff=1, e=0, f=1, gg=0)
Rmodel <- nimbleModel(code, constants, data, inits)


spec <- configureMCMC(Rmodel)

spec$getSamplers()

Rmcmc <- buildMCMC(spec)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
set.seed(0)
Rmcmc$run(3)
Rsamples <- as.matrix(Rmcmc$mvSamples)
set.seed(0)
Cmcmc$run(3)
Csamples <- as.matrix(Cmcmc$mvSamples)

Rnames <- dimnames(Rsamples)[[2]]
Rsamples[, Rnames] - Csamples[, Rnames]
##     x z aaa d dd ddd g
##[1,] 0 0   0 0  0   0 0
##[2,] 0 0   0 0  0   0 0
##[3,] 0 0   0 0  0   0 0
Rsamples[, Rnames]
##             x         z       aaa         d         dd        ddd         g
##[1,] 0.3254630 0.3746455 0.3658422 0.6362147 -0.2698403 -1.1333402 -3.156596
##[2,] 0.1975665 0.3995488 0.3210974 0.3695457 -0.2843177 -0.2239394 -2.711577
##[3,] 0.1900000 0.4199454 0.2194685 0.5430496 -0.9140867  0.1178770 -2.233141


## testing of new dynamic conjugate samplers
### demo2 of check conjugacy
library(nimble)
code <- BUGScode({
    x ~ dbeta(3, 13)
    y[1] ~ dbin(x, 10)
    y[2] ~ dbin(x, 20)
})
data = list(y = c(3,4))
constants <- list()
inits <- list()

Rmodel <- nimbleModel(code, constants, data, inits)

spec <- configureMCMC(Rmodel)
spec$getSamplers()
Rmcmc <- buildMCMC(spec)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Rmcmc$run(10)
Rsamples <- as.matrix(Rmcmc$mvSamples)
set.seed(0)
Cmcmc$run(10)
Csamples <- as.matrix(Cmcmc$mvSamples)

Rsamples - Csamples
as.numeric(Rsamples) - c(0.195510839527966, 0.332847482503424,0.247768152764931, 0.121748195439553, 0.157842271774841, 0.197566496350904, 0.216991517500577, 0.276609942874852, 0.165733872345582, 0.144695512780252)



## testing of new dynamic conjugate samplers
### checkConjugacy_demo3_run.R - various conjugacies
library(nimble)
code <- BUGScode({
    x ~ dgamma(1, 1)       # should satisfy 'gamma' conjugacy class
    a  ~ dnorm(0, x)     # should satisfy 'norm' conjugacy class
    a2 ~ dnorm(0, tau = 3*x+0)
    b  ~ dpois(0+5*x)
    b2 ~ dpois(1*x*1)
    c ~ dgamma(1, 7*x*5)
    for(i in 2:3) {
        jTau[i] <- 1
        jNorm[i] ~ dnorm(c * (a+3) - i, var = jTau[i])
        kTauSd[i] <- 2
        kLogNorm[i] ~ dlnorm(0 - a - 6*i, kTauSd[i])
    }
})
data = list()
constants <- list()
inits <- list()
Rmodel <- nimbleModel(code, constants, data, inits)

spec <- configureMCMC(Rmodel, monitors = c('x', 'c'), control = list(scale = 0.01))
spec$getSamplers()

Rmcmc <- buildMCMC(spec)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Rmcmc$run(10)
Rsamples <- as.matrix(Rmcmc$mvSamples)
set.seed(0)
Cmcmc$run(10)
Csamples <- as.matrix(Cmcmc$mvSamples)

Rnames <- dimnames(Rsamples)[[2]]
Rsamples[, Rnames] - Csamples[, Rnames]

Rsamples[, Rnames] - cbind(c(3.950556165467749, 1.556947815895538, 1.598959152023738, 2.223758981790340, 2.386291653164086, 3.266282048060261, 3.064019155073057, 3.229661999356182, 1.985990552839427, 2.057249437940977), c( 0.010341199485849559, 0.010341199485849559, 0.003846483017887228, 0.003846483017887228, 0.007257679932131476, 0.009680314740728335, 0.012594777095902964, 0.012594777095902964, 0.018179641351556003, 0.018179641351556003))


## testing of new dynamic conjugate samplers
### Dirichlet-multinomial conjugacy
# as of v0.4, exact numerical results here have changed because
# ddirch now sometimes returns NaN rather than -Inf (when an
# alpha is proposed to be negative) -- this changes the RNG
# sequence because NaN values result in no runif() call in decide()
# single multinomial
library(nimble)
set.seed(0)
n <- 100
alpha <- c(10, 30, 15, 60, 1)
K <- length(alpha)
p <- c(.12, .24, .09, .54, .01)
y <- rmulti(1, n, p)
code <- quote({
    y[1:K] ~ dmulti(p[1:K], n);
    p[1:K] ~ ddirch(alpha[1:K]);
    for(i in 1:K) {
        alpha[i] ~ dgamma(.001, .001);
    }
               })
inits <- list(p = rep(1/K, K), alpha = rep(K, K))
constants <- list(n=n, K=K)
data <- list(y = y)
Rmodel <- nimbleModel(code, constants, data, inits)

spec <- configureMCMC(Rmodel, monitors = c('alpha', 'p'))
spec$getSamplers()
Rmcmc <- buildMCMC(spec)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Cmcmc$run(10000)
samples <- as.matrix(Cmcmc$mvSamples)

apply(samples, 2, mean)[6:10]
p
apply(samples, 2, mean)[6:10] - p



## testing of new dynamic conjugate samplers
## rats example
library(nimble)
lst <- readBUGSmodel('rats', dir = getBUGSexampleDir('rats'), returnModelComponentsOnly = TRUE)
code <- lst$model
data <- lst$data[c('Y')]
constants <- lst$data[c('N', 'T', 'x')]
inits <- lst$inits

Rmodel <- nimbleModel(code, constants, data, inits)

spec <- configureMCMC(Rmodel)
spec$getSamplers()

Rmcmc <- buildMCMC(spec)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
set.seed(0)
Cmcmc$run(5000)
samples <- as.matrix(Cmcmc$mvSamples)

a <- samples
b <- samples
a
b
a-b
any(a-b != 0)




## phase transition regression example for Mark.

## generate data
sigma <- 5
a <- 2
b <- 0.5
t0 <- 30
x <- 1:100
n <- length(x)
y <- rnorm(n, a + (x>t0) * (b*(x-t0)), sd=sigma)
plot(x, y)

## fit NIMBLE model
library(nimble)

code <- nimbleCode({
    a ~ dunif(-100, 100)
    b ~ dunif(-100, 100)
    t0 ~ dunif(0, 100)
    sigma ~ dunif(0, 100)
    for(i in 1:n) {
        phase[i] <- step(x[i] - t0)
        mu[i] <- a + phase[i] * b * (x[i]-t0)
        y[i] ~ dnorm(mu[i], sd = sigma)
    }
})
constants <- list(n=n, x=x)
data <- list(y=y)
inits <- list(sigma=1, a=0, b=1, t0=50)

## MCMC Suite will let us do a few MCMCs quickly, and look at results easily
out <- MCMCsuite(
    code, constants, data, inits,
    MCMCs = c('nimble', 'nimble_slice', 'bugs'),
    niter = 20000,
    burnin = 2000,
    makePlot = TRUE,
    savePlot = FALSE,
    calculateEfficiency = TRUE
)

out$summary
out$efficiency
dim(out$samples)
out$timing





## nimble problem with using 'sd' name for 'tau' parameter, e.g.
library(nimble)

code <- nimbleCode({
    sd ~ dunif(0, 1)
    y ~ dnorm(0, tau = sd)
    shape1 ~ dnorm(0, 1)
    f ~ dbeta(mean = y, sd = shape1)
})

m <- nimbleModel(code, inits = list(sd=.1))

m$sd

m$getModelDef()$printDI()

code <- nimbleCode({
    sd2 ~ dunif(0, 1)
    y ~ dnorm(0, tau = sd2)
})
m <- nimbleModel(code, debug=TRUE)

m$sd2






     
## Plot the fake MCMC output
denoverplot(fakemcmc, fakemcmc2)
     denoverplot(fakemcmc, fakemcmc2, style="plain",
                 col=mcmcplotsPalette(3, type="grayscale"),
                 ci=0.95, greek=TRUE)
     denoverplot(fakemcmc, fakemcmc2,
                 plot.title="Comparison of densities of fake data")
     denoverplot(fakemcmc, fakemcmc2,
                 plot.title="Comparison of densities of fake data", greek=TRUE)

setwd('/Users/dturek/GitHub/legacy/autoBlock/')
getwd()

source("autoBlock.R")
load(file.path("data", "model_test.RData"))
saveSamples <- TRUE
niter <- 10000
ab <- autoBlock(code, constants, data, inits, niter, runList, saveSamples = saveSamples)
dftest <- ab$summary
save(dftest, file = file.path("results_samples", "results_test.RData"))
if (saveSamples) {
    burnedSamplesList <- ab$samples
    for (i in 1:length(burnedSamplesList)) burnedSamplesList[[i]] <- burnedSamplesList[[i]][(floor(niter/2) + 1):niter, ]
    save(burnedSamplesList, niter, file = file.path("results_samples", "results_test_samples.RData"))
}

class(burnedSamplesList)
names(burnedSamplesList)

i <- 4
dim(burnedSamplesList[[i]])
dimnames(burnedSamplesList[[i]])

dimnames(ab$samples[[3]])


ret$summary

dftest <- autoBlock(code, constants, data, inits, 10000, runList)$summary
save(dftest, file = file.path("results_hclust_complete2", "results_test.RData"))


rm(list=ls())
getwd()
load('results_samples/results_test.RData')
ls()
dftest

rm(list=ls())
getwd()
load('results_samples/results_test_samples.RData')
ls()
niter
class(burnedSamplesList)
names(burnedSamplesList)
i <- 4
dim(burnedSamplesList[[i]])
dimnames(burnedSamplesList[[i]])







## renaming stupid data frame from the autoBlock results
setwd('~/GitHub/legacy/autoBlock')

rm(list=ls())
file <- 'results/results_samplingEfficiency.RData'
load(file)
ls()
dfsamplingEfficiency <- dfSamplingEfficiency
rm(dfSamplingEfficiency)
ls()
save(list='dfsamplingEfficiency', file = file)

rm(list=ls())
file <- 'results/results_computationalRequirement.RData'
load(file)
ls()
dfcomputationalRequirement <- dfComputationalRequirement
rm(dfComputationalRequirement)
ls()
save(list='dfcomputationalRequirement', file = file)





## working on the spatial capture-recapture (SCR) models

tr<-seq(15,85, length=10)

X<-cbind(rep(tr,each=length(tr)),rep(tr,times=length(tr))) # 100 coord. traps

plot(X, xlim=c(0,100), ylim=c(0,100), pch=3, cex=0.75)

set.seed(10)
xlim <- c(0,100); ylim <- c(0,100)      # Area 100*100=1e4
A <- (xlim[2]-xlim[1])*(ylim[2]-ylim[1])/10000
mu <- 50                 # Density
N <- rpois(1, mu*A) ;N   # Generate population


s <- cbind(runif(N, xlim[1], xlim[2]), runif(N, ylim[1], ylim[2]))
points(s, pch=16, col=2)

sigma <- 5
lambda0 <- 0.4
J <- nrow(X)
K <- 5
yy <- array(NA, c(N, J, K))
for(j in 1:J) {
    dist <- sqrt((X[j,1]-s[,1])^2 + (X[j,2]-s[,2])^2)
    lambda <- lambda0*exp(-dist^2/(2*sigma^2))
    for(k in 1:K) {
        yy[,j,k] <- rpois(N, lambda)
    }
}

n <- apply(yy, c(2,3), sum)

# Plot capture events
tot<-apply(n, 1,sum)
symbols(X, circles=tot, inches=F, bg="#00000022", add=T)
points(X, pch=3, cex=0.75); points(s, pch=16, col=2)


library(nimble)

## define the model
code <- nimbleCode({
    sigma ~ dunif(0,10)
    lam0 ~ dunif(0,5)
    psi ~ dbeta(1,1)
    for(i in 1:M) {
        z[i] ~ dbern(psi)
        s[i,1] ~ dunif(xlim[1], xlim[2])
        s[i,2] ~ dunif(ylim[1], ylim[2])
        for(j in 1:J) {# Number of traps
            dist[i,j] <- (s[i,1] - X[j,1])^2 + (s[i,2] - X[j,2])^2
            lam[i,j] <- lam0*exp(-dist[i,j]/(2*sigma^2))*z[i]
        }
    }
    for(j in 1:J){
        bigLambda[j] <- sum(lam[1:M,j])
        for(k in 1:K) {
            n[j,k] ~ dpois(bigLambda[j])
        }
    }
    N <- sum(z[1:M])
})


M<-200

constants <- list(M = M, K=K, J=J)
n1<-apply(n,1,sum)
data<-list(n=n, X=X, xlim=xlim, ylim=ylim)
s<-cbind(runif(M, xlim[1], xlim[2]), runif(M, ylim[1], ylim[2]))
z<-rep(1,M)
inits <- list (sigma=0.5, lam0=0.1, s=s, z=z)

Rmodel <- nimbleModel(code=code, constants=constants, data=data, inits=inits, check=FALSE) ## check=FALSE is faster

mcmcspec<-configureMCMC(Rmodel, print=TRUE, monitors = c("N", "lam0", "psi", "sigma"))
scrMCMC <- buildMCMC(mcmcspec)
Cmodel <- compileNimble(Rmodel) 
CscrMCMC <- compileNimble(scrMCMC, project = Rmodel)

## It's pretty slow to run so I just want an execution time first for a small sample
t1_100 <- system.time(CscrMCMC$run(100))

## And now I want a decent sample
t1_20k <- system.time(CscrMCMC$run(20000))
t1_20k_samples <- as.matrix(CscrMCMC$mvSamples)
save(t1_100, t1_20k, t1_20k_samples, file = "case1results.Rdata")


## checking dcat distribution for Chris

library(nimble)

code <- nimbleCode({
    y ~ dcat((p[1:2,1:2] %*% ones[1:2])[1:2,1])
})
constants <- list(ones = c(1,1))
data <- list()
inits <- list(y = 2, p = diag(c(.5,.5)))

Rmodel <- nimbleModel(code, constants, data, inits)




library(nimble)
Rmodel <- readBUGSmodel('rats', dir = '~/GitHub/nimble/nimble/packages/nimble/inst/classic-bugs/vol1/rats/')

customSpec <- configureMCMC(Rmodel)
customSpec$getSamplers()
customSpec$removeSamplers(1:65)
customSpec$getSamplers()
customSpec$addSampler(target = 'alpha.c', type = 'slice')
customSpec$addSampler(target = 'beta.c', type = 'slice')
customSpec$addSampler(target = 'tau.c', type = 'slice')
customSpec$addSampler(target = 'tau.alpha', type = 'slice')
customSpec$addSampler(target = 'tau.beta', type = 'slice')
customSpec$getSamplers()
alphaNames <- Rmodel$getDependencies('tau.alpha', self = FALSE, stochOnly = TRUE)
alphaNames
length(alphaNames)  ## 30
for(aN in alphaNames) customSpec$addSampler(target = aN, type = 'slice')
customSpec$getSamplers()
betaNames <- Rmodel$getDependencies('tau.beta', self = FALSE, stochOnly = TRUE)
betaNames
length(betaNames)  ## 30
for(bN in betaNames) customSpec$addSampler(target = bN, type = 'slice')
customSpec$getSamplers()
customSpec


## making a simple MCMC traceplot, and density plot

n <- 1000
x <- rnorm(n, 3, 1)
dev.new(width=4, height=3)
plot(1:n, x, type='l')
dev.copy2pdf(file = '~/Downloads/traceplot.pdf')
dev.new(width=4, height=3)
plot(density(x))
dev.copy2pdf(file = '~/Downloads/densityplot.pdf')


## make sure I get the empirical covariance right

x <- matrix(c(1,2,1,2,1,4,6,5,7,4), 5, 2)
timesRan <- 5
statSums <- apply(x, 2, sum)
statSums <- t(statSums)
statProds <- matrix(0, 2, 2)
for(i in 1:timesRan) statProds <- statProds + t(t(x[i,])) %*% t(x[i,])
statSums
statProds
(statProds - (t(statSums) %*% statSums)/timesRan) / (timesRan-1)
cov(x)
colMeans <- apply(x, 2, mean)
xCentered <- x - rep(colMeans, each=timesRan)
xCentered
(t(xCentered) %*% xCentered) / (timesRan-1)

set.seed(0)
library(nimble)
code <- nimbleCode({
    x[1:d] ~ dmnorm(mu[1:d], cov = Sigma[1:d, 1:d])
})
d <- 4
mu <- 1:d
ch <- array(rnorm(d^2), c(d,d))
Sigma <- ch %*% t(ch)
print(Sigma)
constants <- list(d=d, mu=mu, Sigma=Sigma)
data <- list()
inits <- list(x = rep(1, d))
Rmodel <- nimbleModel(code, constants, data, inits)

spec <- configureMCMC(Rmodel, nodes = NULL)
spec$addSampler('x', 'RW_block', print = FALSE)
spec$addSampler('x', 'RW_block_NEW', print = FALSE)
spec$getSamplers()
Rmcmc <- buildMCMC(spec)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
Cmcmc$run(1000)
samples <- as.matrix(Cmcmc$mvSamples)
cov(samples)
           x[1]       x[2]     x[3]       x[4]
x[1]  2.0934846 -0.1465907 1.107342  1.1252451
x[2] -0.1465907  4.1696671 1.644427 -0.5950559
x[3]  1.1073420  1.6444267 2.581730  1.5097074
x[4]  1.1252451 -0.5950559 1.509707  1.8741393
samples[c(1,250,500,750,1000),]
           x[1]      x[2]       x[3]     x[4]
[1,]  1.0000000  1.000000  1.0000000 1.000000
[2,]  0.2162332  2.260594  1.9804354 2.913190
[3,]  1.8118014  3.966497  2.7563534 2.552825
[4,] -1.6015689 -1.004725 -0.9864974 1.065231
[5,]  0.8037007  2.496745  1.2659404 2.208456






## looking into RStudio hanging problem (once again... Nov 2015)

library(nimble)
code <- nimbleCode({
  a ~ dnorm(0, 1)
  sd ~ dunif(0, 100)
  for(i in 1:N) {
    b[i] ~ dnorm(a, sd = sd)
  }
})
N <- 100
constants <- list(N = N)
data <- list()
inits <- list(a = 0, sd = 1, b = rep(0, N))


Rmodel <- nimbleModel(code, constants, data, inits)


## NEED TO FIND THE PROBLEM IN MODEL CHECKING

## EXAMPLE 1:
rm(list = ls())
library(nimble)
classicDyesModel <- readBUGSmodel('dyes', dir = getBUGSexampleDir('dyes'))

## EXAMPLE 2:
set.seed(0)
library(nimble)
n_obs <- 50
### Parameters
mu <- 2
rho_a <- 1
rho_b <- 0.1
monitored_param <- c('mu', 'rho_a', 'rho_b', 'ngis')
SV =  rgamma(n=n_obs, 1, 1)
### Data generation
lambda <-  exp(mu)*SV
ngis  <-  rpois(n=n_obs, lambda=lambda )
Y <- rgamma(n=n_obs, ngis*rho_a, rho_b)
### Remove zero values from ngis
ngis <- ngis[which(ngis>0)]
abs <- which(Y==0)
n_abs <- length(abs)
pres <- which(Y>0)
n_pres <- length(pres)

# Install from source for OS X, Linux, or Windows:
#install.packages("nimble", repos = "http://r-nimble.org", type = "source")
M1_CPG_nimbleCode <- nimbleCode(
    {
        ## Covariates
        for (s in 1:n_obs){
            lambda[s] <- exp(mu)*SV[s]
        }
        ## Observation model
        ## Strictly positive Observation
        for ( s in 1:n_pres){
            ngis[s] ~ T(dpois(lambda[pres[s]]), 1, )  
            Y[pres[s]] ~ dgamma(rho_a * ngis[s], rho_b)
        }
        ## Zero Observation
        for ( s in 1:n_abs){
            proba[s] <- 1 - exp(-lambda[abs[s]])
            Y[abs[s]] ~ dbern(proba[s])
        }
        ## Prior
        rho_a ~ dgamma(0.01, 0.01)
        rho_b ~ dgamma(0.01, 0.01)
        mu ~ dnorm(0, 0.01)
    })

constants_list <- list(n_obs = length(Y), n_abs=n_abs, n_pres=n_pres, abs=abs, pres=pres)
nimble_data <- list(Y = Y, SV=SV)
## Initial values
inits_fc <- function(chain_id = 1){
  mu <- rnorm(1,0,1)
  rho_a <- runif(1,0,10)
  rho_b <- runif(1,0,10)
  ngis <- round(runif(n_pres, 2, 5), digits=0)
  list(mu=mu, rho_a=rho_a, rho_b=rho_b ,ngis=ngis)
}
init_ll <- lapply(1, function(id) inits_fc(chain_id = id))

Rmodel <- nimbleModel(code= M1_CPG_nimbleCode, name= 'M1_CPG_nimble', constants = constants_list, data =nimble_data, inits = init_ll[[1]], check = TRUE)

calculate(Rmodel, 'ngis[1]')

Rmodel$ngis

check <- TRUE
check <- FALSE

md <- nimbleModel(code= M1_CPG_nimbleCode, name= 'M1_CPG_nimble', constants = constants_list, data =nimble_data, inits = init_ll, returnDef = TRUE)

M1_CPG_nimble <- md$newModel(data =nimble_data, inits = init_ll, check = check)

debug(md$newModel)
md$newModel(data =nimble_data, inits = init_ll, check = TRUE)

debug(model$check)





## random stuff for CR paper....

rm(list=ls())
library(nimble)
##source('~/GitHub/legacy/dipper/dipperCode.R')

Rmodel <- nimbleModel(code_dipper, constants, data, inits)
lnodes <- Rmodel$expandNodeNames('x')
length(lnodes) - nind

rm(list=ls())
trunc <- FALSE
source('~/GitHub/userDistMCMC/defs.R')
source('~/GitHub/userDistMCMC/create_data.R')
Rmodel <- nimbleModel(orchidDHMM$code, orchidDHMM$constants, orchidDHMM$data, orchidDHMM$inits)
topN <- Rmodel$getNodeNames(topOnly = TRUE, stochOnly = TRUE)
topN
length(topN)

con <- orchidDHMM$constants
con$f
con$nind
#### actually 250 individuals in original dataset!!!

con$f
12 - con$f
sum(12 - con$f)

rm(list=ls())
load('~/GitHub/userDistMCMC/results.RData')
or <- results$df[results$df$model == 'orchid', ]
or
or[or$mcmc %in% c('jags', 'nimbleDHMM2'), ]
or1 <- or[or$param == 's[1]', ]
data.frame(mcmc = or1$mcmc, runtime_minutes = or1$timing/60)

rm(list=ls())
load('~/GitHub/userDistMCMC/results.RData')
go <- results$df[results$df$model == 'goose', ]
go
go1 <- go[go$param == 'p[1]', ]
data.frame(mcmc = go1$mcmc, runtime_minutes = go1$timing/60, runtime_hours = go1$timing/60/60)


rm(list=ls())
trunc <- FALSE
source('~/GitHub/userDistMCMC/defs.R')
source('~/GitHub/userDistMCMC/create_data.R')
Rmodel <- nimbleModel(gooseDHMM$code, gooseDHMM$constants, gooseDHMM$data, gooseDHMM$inits)
topN <- Rmodel$getNodeNames(topOnly = TRUE, stochOnly = TRUE)
topN
length(topN)



## testing the MVN conjugacy test from test-mcmc.R
library(nimble)
set.seed(0)
mu0 = 1:3
Q0 = matrix(c(1, .2, .8, .2, 2, 1, .8, 1, 2), nrow = 3)
Q = solve(matrix(c(3, 1.7, .9, 1.7, 2, .6, .9, .6, 1), nrow = 3))
a = c(-2, .5, 1)
B = matrix(rnorm(9), 3)

##### not currently working - see Perry's email of ~ 10/6/14
## code <- nimbleCode({
##   mu[1:3] ~ dmnorm(mu0[1:3], Q0[1:3, 1:3])
##   y[1:3] ~ dmnorm(asCol(a[1:3]) + B[1:3, 1:3] %*% asCol(mu[1:3]), Q[1:3, 1:3])
## })

code <- nimbleCode({
  mu[1:3] ~ dmnorm(mu0[1:3], Q0[1:3, 1:3])
  y_mean[1:3] <- asCol(a[1:3]) + B[1:3, 1:3] %*% asCol(mu[1:3])
  y[1:3] ~ dmnorm(y_mean[1:3], Q[1:3, 1:3])
})

## Simplest version of model w/o 'a' and 'B'
## a = rep(0,3)
## B = diag(rep(1,3))
## code <- nimbleCode({
##   mu[1:3] ~ dmnorm(mu0[1:3], Q0[1:3, 1:3])
##   y[1:3] ~ dmnorm(mu[1:3], Q[1:3, 1:3])
## })


mu <- mu0 + chol(solve(Q0)) %*% rnorm(3)
# make sure y is a vec not a 1-col matrix or get a dimensionality error
y <- c(a + B%*%mu + chol(solve(Q)) %*% rnorm(3))
data = list(mu0 = mu0, Q0 = Q0, Q = Q, a = a, B = B, y = y)

muQtrue = t(B) %*% Q%*%B + Q0
muMeanTrue = c(solve(muQtrue, crossprod(B, Q%*%(y-a)) + Q0%*%mu0))

constants <- list()
inits <- list()

Rmodel <- nimbleModel(code, constants, data, inits)

spec <- configureMCMC(Rmodel)
spec$getSamplers()
Rmcmc <- buildMCMC(spec)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

##debug(Rmcmc$samplerFunctions$contentsList[[1]]$run)
set.seed(0)
Rmcmc$run(5)
Rsamples <- as.matrix(Rmcmc$mvSamples)
Rsamples

set.seed(0)
Cmcmc$run(5)
Csamples <- as.matrix(Cmcmc$mvSamples)
Csamples

Rsamples - Csamples

Cmcmc$run(10000)
samples <- as.matrix(Cmcmc$mvSamples)
apply(samples, 2, mean)
muMeanTrue

apply(samples, 2, mean) - muMeanTrue


source('~/GitHub/nimble/nimble/packages/nimble/inst/tests/test_utils.R')
##test_mcmc
test_mcmc(model = code, name = 'two-level multivariate normal', data = data, seed = 0, numItsC = 10000,
          results = list(mean = list(mu = muMeanTrue),
                           cov = list(mu = solve(muQtrue))),
          resultsTolerance = list(mean = list(mu = rep(.02,3)),
            cov = list(mu = matrix(.01, 3, 3))))



## working the basketball "streaks" problem emailed out by Aaron Strauss
p <- 0.1
ns <- 1
run <- function(p, ns) {
    out <- numeric()
    n <- 0
    while(n < ns) {
        nmade <- 1
        while(rbinom(1, 1, p) == 1) nmade <- nmade + 1
        out <- c(out, nmade)
        n <- length(out)
    }
    return(out)
}

pOfStreaks <- function(p, streaks, niter=1000) {
    pMissing <- missing(p)
    ns <- length(streaks)
    suc <- 0
    for(i in 1:niter) {
        if(pMissing) p <- runif(1)
        thisStreak <- run(p, ns)
        if(all(thisStreak == streaks)) suc <- suc+1
    }
    prob <- suc/niter
    return(prob)
}

pOfStreaks(p=0.2, streaks=4, niter=100000)
0.2^3*(.8)

p <- 0.7
n <- 4
m <- 7
pOfStreaks(p=p, streaks=c(n,m), niter=1000000)
p^(n+m-2)*(1-p)^2

n <- 2
m <- 1
pOfStreaks(streaks=c(n,m), niter=1000000)
1/12

pOfPlessthanPROB <- function(streaks, PROB, niter=1000) {
    ns <- length(streaks)
    suc <- 0
    ranSoFar <- 0
    while(ranSoFar < niter) {
        p <- runif(1)
        thisStreak <- run(p, ns)
        if(all(thisStreak == streaks)) {
            if(p < PROB) suc <- suc+1
            ranSoFar <- ranSoFar+1
        }
    }
    prob <- suc/niter
    return(prob)
}

pOfPlessthanPROB(streaks=c(2,1), PROB=0.5, niter=20000)
11/16

## continusing the "sreaks" basketball problem, now in NIMBLE

library(nimble)

code <- nimbleCode({
    ##p ~ dunif(0, 1)
    p ~ dbeta(1, 1)
    y[1] ~ dnegbin(prob=p, size=1)
    y[2] ~ dnegbin(prob=p, size=1)
})
constants <- list()
data <- list(y = 0:1)  ## the 'streak' history of {2,1} in the original problem
inits <- list(p = 0.5) ## translates into c(1,0) in terms of the negative-binomial distribution

Rmodel <- nimbleModel(code, constants, data, inits)

spec <- configureMCMC(Rmodel)
spec$getSamplers()
Rmcmc <- buildMCMC(spec)
#Rmcmc$samplerFunctions$contentsList[[1]]$run
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
debug(Rmcmc$run)
debug(samplerFunctions[[1]]$run)
Rmcmc$run(1)
set.seed(0)
Cmcmc$run(3)
Rsamples <- as.matrix(Rmcmc$mvSamples)
Csamples <- as.matrix(Cmcmc$mvSamples)
Rsamples
Csamples
Rsamples - Csamples

apply(samples, 2, mean)

mean(samples[,1] > 0.5)
5/16

out <- MCMCsuite(
    code, constants, data, inits,
    MCMCs = c('nimble', 'nimble_noConj'),
    summaryStats = c('mean', 'median', 'sd', 'function(x) mean(x>0.5)'),
    calculateEfficiency = TRUE,
    makePlot = TRUE,
    savePlot = FALSE,
    niter = 100000
)

out$timing
out$summary
11/16


## creating a really small, blank figure, for use in captions in CR-MCMC paper
dev.new(width=3, height=2)
plot(1, 1)
dev.copy2pdf(file='~/GitHub/nimble/nimblePapers/CR-MCMC/blankFigure.pdf')

## testing the problem with naming a variable 'ans' in a nimbleFunction
library(nimble)

rcFunction1 <- nimbleFunction(
    run = function() {
        declare(ans, double(1, 5))
        for(i in 1:5) ans[i] <- i
        print(ans)
    })

rcFunction1()

rcFunction2 <- nimbleFunction(
    run = function() {
        declare(functionAsList, double(1, 5))
        for(i in 1:5) functionAsList[i] <- i
        print(functionAsList)
    })

rcFunction2()

code <- nimbleCode({
    mu ~ dnorm(0, 1)
    mu2 <- timesTwo(mu)
})

Rmodel <- nimbleModel(code, inits = list(mu=1), check=FALSE)

##undebug(calculate)
##debug(nimble:::rCalcNodes)
##Rmodel$nodes$mu2$calculate
##Rmodel$nodes$mu2$calculate()
##Rmodel$nodes$mu2$simulate
##Rmodel$nodes$mu2$simulate()

calculate(Rmodel)
Cmodel <- compileNimble(Rmodel)

Rmodel$mu
Rmodel$mu2
Cmodel$mu
Cmodel$mu2

Rmodel$mu <- 2
calculate(Rmodel)
Rmodel$mu2


## finding the posterior pairwise correlations of
## the blocked parameters for the orchid model
## (multistate CR-MCMC, userDistMCMC)

library(nimble)
rm(list=ls())
source('~/GitHub/userDistMCMC/defs.R')
load('~/GitHub/userDistMCMC/models.RData')

Rmodel <- nimbleModel(orchidDHMM$code, orchidDHMM$constants, orchidDHMM$data, orchidDHMM$inits)
spec <- configureMCMC(Rmodel)
spec$getSamplers()
spec$removeSamplers(numeric(0))
spec$removeSamplers(c('b[1:2]','c[1:2]'))
spec$addSampler(c('b[1]','b[2]'), 'RW_block')
spec$addSampler(c('c[1]','c[2]'), 'RW_block')
spec$getSamplers()
Rmcmc <- buildMCMC(spec)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
## n <- 10000
n <- 50000
system.time(Cmcmc$run(n)) / 60
samples <- as.matrix(Cmcmc$mvSamples)
samples2 <- samples[(n/2+1):n, ]   ## truncate to second half of samples
apply(samples2, 2, mean)
Cov <- cov(samples2)
Cov
Cor <- cov2cor(Cov)
Cor
sort(as.numeric(abs(Cor)))
Cor['b[1]', 'b[2]']
## [1] 0.9866091
Cor['c[1]', 'c[2]']
## [1] 0.9675666


## printing the numeric performance results for CR-MCMC paper, or userDistMCMC paper
load('~/GitHub/userDistMCMC/results.RData')
##load('~/GitHub/userDistMCMC/resultsNew.RData')
source('~/GitHub/userDistMCMC/defs.R')
options(digits = 3)
df <- results$df
for(mod in unique(df$model)) {
    dfmod <- df[df$model==mod, ]
    message('******************************************')
    message(mod, ' model')
    message('******************************************')
    for(mc in unique(dfmod$mcmc)) {
        dfmodmcmc <- dfmod[dfmod$mcmc==mc, ]
        message(mc, ' MCMC:')
        message('minimum: ',  min(dfmodmcmc$Efficiency))
        message('mean: ',    mean(dfmodmcmc$Efficiency))
    }
    message('******************************************')
}

df[df$model=='orchid',]
df[df$model=='goose',]

results$check()
results$quickplot()


## problem with multivariate conjugate samplers... ?
## fixed.  problem was in vector demoting in the conjugacy definition
library(nimble)
code <- nimbleCode({
    x[1:d] ~ dmnorm(mu[1:d], cov = Sigma[1:d, 1:d])
    y[1:d] ~ dmnorm(x[1:d], cov = Sigma[1:d, 1:d])
})
d <- 2
mu <- rep(0, d)
Sigma <- diag(d)
constants <- list(d=d, mu=mu, Sigma=Sigma)
Rmodel <- nimbleModel(code, constants, data=list(y=rep(0,d)), inits=list(x=rep(0, d)))
Cmodel <- compileNimble(Rmodel)

spec <- configureMCMC(Rmodel)
spec$getSamplers()
Rmcmc <- buildMCMC(spec)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Cmcmc$run(10000)
samples <- as.matrix(Cmcmc$mvSamples)
apply(samples, 2, mean)
solve(cov(samples))

library(nimble)
nimble:::conjugateSamplerDefinitions$sampler_conjugate_dnorm
nimble:::conjugateSamplerDefinitions$sampler_conjugate_dbeta
nimble:::conjugateSamplerDefinitions$sampler_conjugate_dmnorm



## testing different priors for rho in package gpmanagement
library(gpmanagement)
library(ggplot2)
example_data <- function(p, x0, Tobs){
    f <- function (x, h){
        sapply(x, function(x) {
            x <- pmax(0, x - h)
            x * exp(p[1] * (1 - x/p[2]) * (x - p[3])/p[2])
        })
    }
    sigma_g <- 0.1
    pdfn <- function(x, mu, sigma = sigma_g){
        dlnorm(x, log(mu), sdlog = sigma)
    }
    z_g <- function() rlnorm(1, 0, sigma_g)
    x <- numeric(Tobs)
    x[1] <- x0
    for(t in 1:(Tobs-1))
        x[t+1] = z_g() * f(x[t], h=0)
    obs <- data.frame(x = c(0, 
                          pmax(rep(0,Tobs-1), x[1:(Tobs-1)])), 
                      y = c(0, 
                          x[2:Tobs]))
    obs
}
myfit <- function(xObs, yObs, priors = NULL, niter = 100000){    ## default value for priors
  xPred <- seq(0, 1.1 * max(xObs), length = 50)
  fit <- gp_setup(xObs, yObs, xPred, priors)     ## included priors argument
  Cmcmc <- fit$Cmcmc 
  Cpred <- fit$Cpred
  Cmodel <- fit$Cmodel
  system.time(Cmcmc$run(niter))
  samples <- as.matrix(Cmcmc$mvSamples)
  system.time(Cpred$run(samples))
  E <- Cpred$getE()
  V <- sqrt(diag(Cpred$getC()))
  list(samples = tidyr::gather(as.data.frame(samples)),
       pred = data.frame(x = xPred, y = E, ymin = E - V, ymax = E + V))
}

densPlot <- function(fit, paramName, paramPrior) {
        dev.new()
        samp <- fit$samples[fit$samples$key==paramName, 'value']
        plot(density(samp), main=paste0('black posterior, red prior: ', paramName, ' ~ ', deparse(paramPrior)))
        xs <- seq(range(samp)[1], range(samp)[2], length.out = 500)
        yCall <- as.call(c(list(as.name(paramPrior[[1]])), list(quote(xs)), as.list(paramPrior[-1])))
        ys <- eval(yCall, envir = environment())
        lines(xs, ys, col='red')
}

## Carrying capcacity at 10, critical point at 5, x0=6
## 100, 50, x0=60
plotPred <- function(ccap, crit, x0, Tobs=40, plotData=FALSE, rhoPrior, plotRhoDensity=TRUE, niter=100000, plotSigGPDensity=FALSE, plotSigOEDensity=FALSE, plotPredictionIntervals=TRUE) {
    set.seed(0)
    obs <- example_data(p=c(2,ccap,crit), x0=x0, Tobs=Tobs)
    if(plotData) {dev.new(); qplot(seq_along(obs$x[-1]), obs$x[-1]) + geom_line()}
    priors <- expression({
        rho ~ X
        sigGP ~ dunif(0, 1e5)
        sigOE ~ dunif(0, 1e5)
    })
    if(missing(rhoPrior)) stop('must specify a prior for rho')
    rhoPrior[-1] <- lapply(rhoPrior[-1], function(param)
        eval(eval(substitute(substitute(PARAM, list(xObs=obs$x[-1])), list(PARAM=param)))))
    print(rhoPrior)
    priors[[1]][[2]][[3]] <- rhoPrior
    fit <- myfit(obs$x, obs$y, priors = priors, niter=niter)
    if(plotRhoDensity) densPlot(fit, 'rho', rhoPrior)
    if(plotSigGPDensity) densPlot(fit, 'sigGP', priors[[1]][[3]][[3]])
    if(plotSigOEDensity) densPlot(fit, 'sigOE', priors[[1]][[4]][[3]])
    if(plotPredictionIntervals) {
        dev.new()
        print(ggplot2::ggplot(fit$pred) + 
                  geom_ribbon(aes(x = x, y = y, ymin = ymin, ymax = ymax), fill = "grey80") +
                      geom_line(aes(x = x, y = y), size = 1) + 
                          geom_point(data = obs, aes(x, y)) +
                              coord_cartesian(xlim = range(c(obs$x, fit$pred$x)), 
                                              ylim = range(c(obs$y, fit$pred$y))) +
                                                  ggtitle(paste0('rho ~ ', deparse(rhoPrior))))
    }
}

rp <- quote(dunif(0, diff(range(xObs))/sqrt(6)/2))

plotPred(rhoPrior = rp, ccap= 10, crit= 5, x0= 6)
plotPred(rhoPrior = rp, ccap=100, crit=50, x0=60)

plotPred(rhoPrior = rp, ccap=100, crit=50, x0=60, plotSigGPDensity=TRUE, plotSigOEDensity=TRUE)

plotPred(rhoPrior = quote(dunif(0,200)), ccap=10, crit=5, x0=6)

for(uMax in 5:6) {
    rp <- substitute(dunif(0, UMAX), list(UMAX = as.numeric(uMax)))
    plotPred(rhoPrior = rp, ccap=10, crit=5, x0=6, niter=10000)
}



## testing ggplot() inside a loop??
library(ggplot2)
df <- data.frame(a=1:2, b=3:4)
for(i in 1:3) {
    dev.new()
    ## need explict print(...) to make ggplot() plot, when inside a loop!!!
    print(ggplot(df, aes(a, b)) + geom_line() + ggtitle(i))
}



## testing use of substitute...
## which I used to be more comfortable with!

f <- function(a) {
    b <- substitute(a)
    print(b)
    print(class(b))
    vv <- eval(b)
    print(vv)
    print(class(vv))
}

f(dnorm(0,1))
f(4+5)

x <- f(quote(dnorm(0,1)))
x
class(x)


## testing NIMBLE compiler to handle pi
nfDef <- nimbleFunction(
    setup = function() {},
    run = function() {
        a <- pi
        print(a)
        b <- 2 * pi
        print(b)
    }
)
Rnf <- nfDef()
Rnf$run()
Cnf <- compileNimble(Rnf)    ## ERROR
Cnf$run()


## example of trying to use NIMBLE functions in another package
## this uses my minimalist nimtest package; same thing I sent to Duncan
library(devtools)
install_github('danielturek/nimtest')
library(nimtest)
library(nimble)
## First Test:
## a nimbleFunction defined as a package function
nfDefinition
Rnf <- nfDefinition()
Cnf <- compileNimble(Rnf)
Rnf$run()
Cnf$run()
## Second Test:
## A package function which internally defines and returns a nimbleFunction definition
## This one gives ugly looking warnings
nextFunction
def <- nextFunction()        ## WARNINGS
Rnf <- def()
Cnf <- compileNimble(Rnf)    ## WARNINGS
Rnf$run()
Cnf$run()

## alternate.  from shell:
## R CMD BUILD nimtest
## R CMD install nimtest_0.0.0.9000.tar.gz
library(nimtest)
library(nimble)
nfDefinition    ## finds it -- it created the NIMBLE function at build time
Rnf <- nfDefinition()  ## can't use it
Cnf <- compileNimble(Rnf)
Rnf$run()
Cnf$run()
nextFunction
nextFunction()
def <- nextFunction()
Rnf <- def()  ## can't use it
Cnf <- compileNimble(Rnf)
Rnf$run()
Cnf$run()


## testing out the new elliptical slice sampler (ess)
library(nimble)
code <- nimbleCode({
    x[1:d] ~ dmnorm(mu[1:d], cov = Sigma[1:d, 1:d])
})
d <- 2
##mu <- rep(0, d)
mu <- 1:d
##Sigma <- diag(d);  print(Sigma)
ch <- array(c(1, 0.7, 0, 2), c(2,2));  Sigma <- ch %*% t(ch);  print(Sigma)
constants <- list(d=d, mu=mu, Sigma=Sigma)
data <- list()
inits <- list(x = rep(1, d))
Rmodel <- nimbleModel(code, constants, data, inits)
Cmodel <- compileNimble(Rmodel)
spec <- configureMCMC(Rmodel, nodes = NULL)
spec$addSampler('x', 'ess', print=FALSE)
spec$getSamplers()
Rmcmc <- buildMCMC(spec)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0); Cmcmc$run(20000)
samples <- as.matrix(Cmcmc$mvSamples)
apply(samples, 2, mean)
cov(samples)   ## empirical cov seems WRONG!  about half of what it should be?
Sigma


## conjugate MVN-MVN
set.seed(1)
library(nimble)
code <- nimbleCode({
    x[1:d] ~ dmnorm(mu_x[1:d], prec = prec_x[1:d, 1:d])
    y[1:d] ~ dmnorm(x[1:d], prec = prec_y[1:d, 1:d])
})
d <- 3
mu_x <- rnorm(d)
temp <- array(rnorm(d^2), c(d,d))
prec_x <- solve(temp %*% t(temp))
temp <- array(rnorm(d^2), c(d,d))
prec_y <- solve(temp %*% t(temp))
y <- rnorm(d)
constants <- list(d = d, mu_x = mu_x, prec_x = prec_x, prec_y = prec_y)
data <- list(y = y)
inits <- list(x = rep(0, d))
Rmodel <- nimbleModel(code, constants, data, inits)
Cmodel <- compileNimble(Rmodel)

## block sampling
spec <- configureMCMC(Rmodel, nodes = NULL)
spec$addSampler('x', 'RW_block', print = FALSE)
spec$getSamplers()
Rmcmc <- buildMCMC(spec)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
set.seed(0); Cmcmc$run(100000)
samples <- as.matrix(Cmcmc$mvSamples)
solve(prec_x + prec_y)
cov(samples)   ## empirical cov seems WRONG!  about half of what it should be?
solve(prec_x + prec_y) %*% (prec_y %*% y + prec_x %*% mu_x)
apply(samples, 2, mean)

## ess sampling
spec <- configureMCMC(Rmodel, nodes = NULL)
spec$addSampler('x', 'ess', print = FALSE)
spec$getSamplers()
Rmcmc <- buildMCMC(spec)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
set.seed(0); Cmcmc$run(100000)
samples <- as.matrix(Cmcmc$mvSamples)
solve(prec_x + prec_y)
cov(samples)   ## empirical cov seems WRONG!  about half of what it should be?
solve(prec_x + prec_y) %*% (prec_y %*% y + prec_x %*% mu_x)
apply(samples, 2, mean)



## testing elliptical slice sampler on 'spatial' model
library(nimble)
load('~/GitHub/legacy/autoBlock/data/model_spatial.RData')

out <- MCMCsuite(
    code, constants, data, inits,
    MCMCs = c('block', 'ess'),
    MCMCdefs = list(
        block = quote({
            spec <- configureMCMC(Rmodel, nodes = NULL)
            spec$addSampler('mu', 'RW', print=FALSE)
            spec$addSampler(c('sigma', 'rho'), 'RW_block', print=FALSE)
            spec$addSampler('g[1:148]', 'RW_block', print=FALSE)
            spec
        }),
        ess = quote({
            spec <- configureMCMC(Rmodel, nodes = NULL)
            spec$addSampler('mu', 'RW', print=FALSE)
            spec$addSampler(c('sigma', 'rho'), 'RW_block', print=FALSE)
            spec$addSampler('g[1:148]', 'ess', print=FALSE)
            spec
        })),
    niter = 10000,
    burnin = 2000,
    monitors = c('mu', 'sigma', 'rho', 'g'),
    makePlot = FALSE
)

out$timing
out$summary[,,1:5]    ## WRONG again...
cov(t(out$samples['block', -(1:3), ])) - cov(t(out$samples['ess', -(1:3), ]))






## Taken from Carl's lab notebook:
## Gpdp Via MDPtoolbox Cont; 30 Jul 2015
## http://www.carlboettiger.info/2015/07/30/gpdp-via-mdptoolbox-cont.html
## (Markov Decision Processes (MDP) toolbox)
##knitr::opts_chunk$set(eval=FALSE)
##devtools::install_github("cboettig/gpmanagement@b3b765cbceb51c9b0b8cb2724e395353ec365df9")
library("MDPtoolbox")
library("gpmanagement")
library("tidyr")
library("dplyr")
library("ggplot2")

## True model
p <- c(2, 100, 50)
f <- function (x, h){
    sapply(x, function(x) {
        x <- pmax(0, x - h)
        x * exp(p[1] * (1 - x/p[2]) * (x - p[3])/p[2])      #### Allen model
    })
}

sigma_g <- 0.1
pdfn <- function(x, mu, sigma = sigma_g){
    dlnorm(x, log(mu), sdlog = sigma)
}

z_g <- function() rlnorm(1, 0, sigma_g)
set.seed(0)
Tobs <- 40

x <- numeric(Tobs)
x[1] <- 60
for(t in 1:(Tobs-1)) x[t+1] = z_g() * f(x[t], h=0)
obs <- data.frame(x = c(0, pmax(rep(0,Tobs-1), x[1:(Tobs-1)])), 
                  y = c(0, x[2:Tobs]))
xObs <- obs$x
yObs <- obs$y
xPred <- seq(0, 1.1 * max(xObs), length = 50)
qplot(seq_along(x), x) + geom_line()

## Now the GP estimation from NIMBLE. Let’s emphasize shorter length-scales with the prior to compare:

## having this work with output from `nimbleCode` might be more natural than `expression`
priors <- expression({
  rho ~ dgamma(1, 1)
  sigGP ~ dunif(0, 1e5)
  sigOE ~ dunif(0, 1e5)
})

fit <- gp_setup(xObs, yObs, xPred)

Cmcmc <- fit$Cmcmc 
Cpred <- fit$Cpred
Cmodel <- fit$Cmodel
system.time(Cmcmc$run(100000))
samples <- as.matrix(Cmcmc$mvSamples)
## basic sanity check
testthat::expect_identical(Cmodel$getNodeNames(topOnly = TRUE), colnames(samples))

## predict from GP model using posterior MCMC samples
system.time(Cpred$run(samples))

## Posteriors
samples <- as.data.frame(as.matrix(Cmcmc$mvSamples))
df <- tidyr::gather(samples)

ggplot(df) + 
  geom_density(aes(value)) + 
  facet_wrap(~key, scale='free')

## extract predictions: E and C
E <- Cpred$getE()
C <- Cpred$getC()

obs <- data.frame(x = xObs, y = yObs)
pred <- data.frame(x = xPred, y = E, ymin = E - sqrt(diag(C)), ymax = E + sqrt(diag(C)))

ggplot2::ggplot(pred) + 
  geom_ribbon(aes(x = x,y = y, ymin = ymin, ymax = ymax), fill = "grey80") +
  geom_line(aes(x = x, y = y), size=1) + 
  geom_point(data = obs, aes(x,y)) +
  coord_cartesian(xlim = range(c(xObs, xPred)), ylim = range(c(yObs,E))) +
  theme_bw()

## Decision theory
states <- xPred # Vector of all possible states
actions <- states # Vector of actions: harvest
## Let’s consider a slight variation of the most trivial utility function:
## one which explicitly adds a cost to completely exhausting the stock
## (or reducing the stock by more than, say 95% in this case.)
## This should be somewhat similar to the impact of no discount rate.

## Utility function
discount = 0.99

#get_utility <- function(x,h) pmin(x,h)
#R <- outer(states, actions, get_utility)

R <- sapply(actions, function(h){
    sapply(states, function(x){
        if(h < x) h else - 1 * max(states)
    })
})

## Implementing policy
z <- function() rlnorm(1, meanlog = 0, sdlog = sigma_g)

simulate_policy <- function(states, actions, policy, f, z, s0, steps = 50, utility = function(s,a) NA, discount = 1){
    s <- numeric(steps)
    a <- numeric(steps)
    u <- numeric(steps)
    s[1] <- s0
    for(t in 1:(steps-1)){
        a[t] <- actions[policy[which.min(abs(states - s[t]))]]
        s[t+1] <- z() * f(s[t], a[t])
        u[t] <- utility(s[t], a[t]) * discount ^ t
    }
    ## Final action determined but not implemented
    a[steps] <- actions[policy[which.min(abs(states - s[t]))]]
    data.frame(time = 1:steps, state = s, action = a, utility = u)
}

## GP model
gp_matrix <- function(states, actions, E, C){
    transition <- array(0, dim = c(length(states), length(states), length(actions)))
    K <- length(states)
    sigmas <- sqrt(diag(C))
    for (k in 1:length(states)) {
        for (i in 1:length(actions)) {
            nextpop <- E[k] - actions[i]
            if(nextpop <= 0) {
                transition[k, , i] <- c(1, rep(0, K - 1))
            } else {
                transition[k, , i] <- dnorm(states, nextpop, sigmas[i]) / sum(dnorm(states, nextpop, sigmas[i]))
            }
        }
    }
    transition
}

P_gp <- gp_matrix(states, actions, E, C)
mdp_check(P = P_gp, R = R)

gp <- mdp_value_iteration(P_gp, R, discount = discount, epsilon = 0.00001, max_iter = 5e3, V0 = numeric(length(states)))

plot(states, states - actions[gp$policy],  xlab="Population size", ylab="Escapement")

data.frame(reps = 1:50) %>% 
    group_by(reps) %>% 
        do(simulate_policy(states, actions, gp$policy, f, z, s0 = 100, steps = 20, utility = pmin, discount = discount)[c("time", "state", "utility")]) -> sims 

mean(sims$utility)
ggplot(sims) + geom_line(aes(time, state, group = reps), alpha = 0.3, col = "darkblue")
## With this amount of data, the gp solution is too cautious, and avoids any exploitation.

## Simulate under the true model

data.frame(reps = 1:50) %>% 
    group_by(reps) %>% 
        do(simulate_policy(states, actions, gp$policy, f, z, s0 = 100, steps = 20, utility = pmin, discount = discount)[c("time", "state", "utility")]) -> sims 

mean(sims$utility)
## (Average utility is approximate here since it does not include penalty;
## since a function and not a matrix is requred by this function at this time.)

ggplot(sims) + geom_line(aes(time, state, group = reps), alpha = 0.3, col = "darkblue")

P <- transition_matrix(states, actions, f, pdfn)
mdp_check(P = P, R = R)

mdp <- mdp_value_iteration(P, R, discount = discount, epsilon = 0.001, max_iter = 5e3, V0 = numeric(length(states)))
plot(states, states - actions[mdp$policy],  xlab="Population size", ylab="Escapement")

## Note that the altered award structure has almost no effect on the optimal policy
## given the true model, other than to avoid harvesting directly to zero even when
## the stock cannot persist, due to the explicit penalty for doing so.

data.frame(reps = 1:50) %>% 
    group_by(reps) %>% 
        do(simulate_policy(states, actions, mdp$policy, f, z, s0 = 100, steps = 20, utility = pmin, discount = discount)[c("time", "state", "utility")]) -> sims 

mean(sims$utility)

ggplot(sims) + geom_line(aes(time, state, group = reps), alpha = 0.3, col = "darkblue")





## Taken from Carl's lab notebook:
## MDPtoolbox Allen Model; 29 Jul 2015
## http://www.carlboettiger.info/2015/07/29/mdptoolbox-allen-model.html
## (Markov Decision Processes (MDP) toolbox)
library("MDPtoolbox", quietly = TRUE)
library("ggplot2", quietly = TRUE)
K <- 10 # state space limit
states <- 0:K # Vector of all possible states
actions <- states # Vector of actions: harvest

sigma_g = 0.1
p <- c(2, 15, 5)

f <- function (x, h){
    sapply(x, function(x) {
        x <- pmax(0, x - h)
        x * exp(p[1] * (1 - x/p[2]) * (x - p[3])/p[2])      #### Allen model
    })
}

pdfn <- function(x, mu, sigma = sigma_g){
    dlnorm(x, log(mu), sdlog = sigma)
}

## Utility function
discount = 0.95
get_utility <- function(x,h) {
    pmin(x,h)
}

R <- outer(states, actions, get_utility)

transition_matrix <- function(states, actions, f, pdfn){
    ## Initialize
    transition <- array(0, dim = c(length(states), length(states), length(actions)))
    K <- length(states)
    for (k in 1:length(states)) {
        for (i in 1:length(actions)) {
            ## Calculate the transition state at the next step, given the 
            ## current state k and action i (harvest H[i])
            nextpop <- f(states[k], actions[i])
            ## Population always extinct if this is negative.
            ## since multiplicitive shock z_t * f(n) < 0 for all f(n) < 0
            if(nextpop <= 0)
                transition[k, , i] <- c(1, rep(0, length(states) - 1))
            ## Implement demographic stochasticity 
            else {
                ## Cts distributions need long-tailed denominator as normalizing factor:
                fine_states <- seq(min(states), 10 * max(states), by = states[2] - states[1])
                N <- sum(pdfn(fine_states, nextpop))  
                transition[k, , i] <-pdfn(states, nextpop) / N
                ## We need to correct this density for the final capping state ("Pile on boundary")
                ## (discrete or cts case)
                ## this can be a tiny but negative value due to floating-point errors.
                ## so we take max(v,0) to avoid
                transition[k, K, i] <- max(1 - sum(transition[k, -K, i]), 0)
            }
        } 
    }
    transition
}

P <- transition_matrix(states, actions, f, pdfn)
apply(P, c(1,3), sum)  ## double-check

## Using toolbox
mdp_check(P = P, R = R)
mdp <- mdp_value_iteration(P, R, discount = discount, epsilon = 0.001, max_iter = 5e3, V0 = numeric(length(states)))
plot(states, states - actions[mdp$policy],  xlab="Population size", ylab="Escapement")

## Compare to Reed
## From Reed (1979) we know that the optimal solution is a constant-escapement rule
## when the growth function in convex. Note that this condition is violated by the
## growth function with alternative stable states (Allen/Ricker-Allee model),
## resulting in a very different optimal policy:
## f^prime(s^star) = 1/discount
## For growth-rate function f, where discount is the discount factor and s∗ the stock size
## for the constant escapement. Analytic solutions are clearly possible for certain
## growth functions, but here I’ve just implemented a generic numerical solution.
fun <- function(x) - f(x,0) + x / discount
out <- optimize(f = fun, interval = c(0,K))
S_star <- out$minimum
exact_policy <- sapply(states, function(x) if(x < S_star) 0 else x - S_star)
D <- actions[mdp$policy]
# The difference between Bellman and the analytical solution is small:
plot(states, states - D,  xlab="Population size", ylab="Escapement", ylim = c(0, 1.2*max(states-D)))
lines(states, states - exact_policy)




## Carl's examples using MDPtoolbox
## MDPtoolbox Ex 2; 14 Jul 2015
## (Markov Decision Processes (MDP) toolbox)
## http://www.carlboettiger.info/2015/07/14/mdptoolbox-ex-2.html
## Adapted from Marescot et al. appendix 5, to Reed optimal control problem, including direct comparison against (semi) analytic optimum.

## step 1: define objectives
## This is a conceptual step which does not require coding
## step 2: define states
##K <- 150 # state space limit
K <- 10 # state space limit
states <- 0:K # Vector of all possible states

## step 3: define control actions
## Vector of actions: harvest
H <- states

## step 4: define dynamic model (with demographic parameters)
##p <- c(6,0.05)
p <- c(3, 0.2)
f <- function(x, h){
  A <- p[1] 
  B <- p[2] 
  s <- pmax(x-h, 0)
  A * s/(1 + B * s)   #### Bellman equation
}
sigma_g = 0.1

## step 5: define utility
## Utility function
get_utility <- function(x,h) {
    pmin(x,h)
}

## step 6: solve bellman equation with value iteration
## Initialize transition matrix
transition <- array(0, dim = c(length(states), length(states), length(H)))
## Initialize utility matrix
utility <- array(0, dim = c(length(states), length(H)))
## Fill in the transition and utility matrix
## Loop on all states
for (k in 0:K) {
    ## Loop on all actions
    for (i in 1:length(H)) {
        ## Calculate the transition state at the next step, given the 
        ## current state k and the harvest H[i]
        nextpop <- f(k, H[i])
        if(nextpop <= 0)
            transition[k+1, , i] <- c(1, rep(0, length(states) - 1))
        ## Implement demographic stochasticity by drawing probability from a density function
        else {
            ## We need to correct this density for the final capping state ("Pile on boundary").
            ## For discrete probability distribution,
            ## this is easy if `states` includes all possible
            ## discrete states below the capping state.
            ## (e.g. all non-negative integers less than K).  
            ## For a continuous distribution, this is more problematic
            ## as we have to first normalize the densities.
            ## EDIT: this can be negative, due to floating-point errors.
            ## so we take max(v,0) to avoid.
            ## Get long-tailed denominator as normalizing factor (continuous distributions only):
            fine_states <- seq(min(states), 10 * max(states), by = states[2]-states[1])
            N <- sum(dlnorm(fine_states, log(nextpop), sdlog = sigma_g))
            transition[k+1, , i] <- dlnorm(states, log(nextpop), sdlog = sigma_g) / N
            ## We need to correct this density for the final capping state ("Pile on boundary")
            transition[k+1, K+1, i] <- max(1 - sum(transition[k+1, -(K+1), i]), 0)
        }
        ## Compute utility
        utility[k+1, i] <- get_utility(k, H[i])
    } # end of action loop
} # end of state loop

## Solution calculated explicitly:
## The backward iteration consists in storing action values in the vector
## Vt which is the maximum of utility plus the future action values for all
## possible next states. Knowing the final action values, we can then backwardly
## reset the next action value Vtplus to the new value Vt.  We start The backward
## iteration at time T-1 since we already defined the action value at Tmax.

## Discount factor
discount <- 0.95
## Action value vector at tmax
Vtmax <- numeric(length(states))
## Action value vector at t and t+1
Vt <- numeric(length(states))
Vtplus <- numeric(length(states))
## Optimal policy vector
D <- numeric(length(states))
## Time horizon
##Tmax <- 150
Tmax <- 4

for (t in (Tmax - 1):1) {
    ## We define a matrix Q that stores the updated action values for 
    ## all states (rows)
    ## actions (columns)
    Q <- array(0, dim = c(length(states), length(H)))
    for (i in 1:length(H)) {
        ## For each harvest rate we fill for all states values (row) 
        ## the ith column (Action) of matrix Q
        ## The utility of the ith action recorded for all states is 
        ## added to the product of the transition matrix of the ith 
        ## action by the action value of all states 
        Q[,i] <- utility[, i] + discount * (transition[,,i] %*% Vtplus)
    } # end of the harvest loop
    ## Find the optimal action value at time t is the maximum of Q
    Vt <- apply(Q, 1, max)
    ## After filling vector Vt of the action values at all states, we 
    ## update the vector Vt+1 to Vt and we go to the next step standing 
    ## for previous time t-1, since we iterate backward
    Vtplus <- Vt
} # end of the time loop

## Find optimal action for each state
for (k in 0:K) {
    ## We look for each state which column of Q corresponds to the 
    ## maximum of the last updated value 
    ## of Vt (the one at time t + 1). If the index vector is longer than 1 
    ## (if there is more than one optimal value we chose the minimum 
    ## harvest rate)
    D[k + 1] <- H[(min(which(Q[k + 1, ] == Vt[k + 1])))]
}

## plot solution
plot(states, states - D, xlab="Population size", ylab="Escapement")

## proof of optimality: compare with analytical solution
fun <- function(x) - f(x,0) + x / discount
out <- optimize(f = fun, interval = c(0,K))
S_star <- out$minimum
exact_policy <- sapply(states, function(x) if(x < S_star) 0 else x - S_star)

## The difference between Bellman equation solution and the analytical solution is small:
plot(states, states - D, xlab="Population size", ylab="Escapement", ylim = c(0, 1.2*max(states-D)))
lines(states, states - exact_policy)

## Using MDPtoolbox
library('MDPtoolbox')
mdp_check(P = transition, R = utility)
out <- mdp_value_iteration(transition, utility, discount = discount, epsilon = 0.001, max_iter = 5e3, V0 = Vtmax)
plot(states, states - D, xlab="Population size", ylab="Escapement")
lines(states, states - H[out$policy], col="red", lty=2)




## testing using new NIMBLE option: verbose
nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)
nimbleOptions('verbose')

library(nimble)

code <- nimbleCode({
    a ~ dnorm(0, 1)
})
constants <- list()
data <- list()
inits <- list(a = 0)

Rmodel <- nimbleModel(code, constants, data, inits)


## testing the new includeEfficiency option for MCMCsuite
## calculates N, ESS, and Eff

library(nimble)
load('~/GitHub/legacy/autoBlock/data/model_litters.RData')
load('~/GitHub/legacy/autoBlock/data/model_ice.RData')

out <- MCMCsuite(code, constants, data, inits,
##                 MCMCs = c('nimble', 'nimble_slice', 'autoBlock'),
                 MCMCs = c('nimble', 'nimble_slice'),
                 makePlot = FALSE,
                 calculateEfficiency = TRUE
                 )

out$timing

out$summary

apply(out$summary[, 'efficiency', ], 1, min)
apply(out$summary[, 'efficiency', ], 1, mean)

out$efficiency


## running OpenBUGS installation from OSX (using wine)

##library(R2OpenBUGS)
library(R2WinBUGS)
library(BRugs)

## schools data in the R2OpenBUGS library
data(schools)
schools

## define the model
nummodel <- function() {
    for (j in 1:J) {
        y[j] ~ dnorm(theta[j], tau.y[j])
        theta[j] ~ dnorm(mu.theta, tau.theta)
        tau.y[j] <- pow(sigma.y[j], -2)
    }
    mu.theta ~ dnorm(0, 1E-6)
    tau.theta <- pow(sigma.theta, -2)
    sigma.theta ~ dunif(0, 1000)
}

## write the model code out to a file
write.model(nummodel, 'nummodel.txt')
model.file1 = paste(getwd(), 'nummodel.txt', sep='/')
## and let's take a look:
file.show('nummodel.txt')

## prepare the data for input into OpenBUGS
J <- nrow(schools)
y <- schools$estimate
sigma.y <- schools$sd
data <- list(J = J, y = y, sigma.y = sigma.y)

## initialization of variables
inits <- function() {
    list(theta = rnorm(J, 0, 100), mu.theta = rnorm(1, 0, 100), sigma.theta = runif(1, 0, 100))
}

## set the WINE working directory and the directory to OpenBUGS
## change the OpenBUGS.exe location as necessary
WINE <- '/opt/local/bin/wine'
WINEPATH <- '/opt/local/bin/winepath'
##OpenBUGS.pgm <- '/Users/dturek/.wine/drive_c/Program Files/OpenBUGS/OpenBUGS323/OpenBUGS.exe'
bugs.directory <- '/Users/dturek/.wine/drive_c/Program Files/OpenBUGS/OpenBUGS323'

## these are the parameters to save
parameters = c('theta', 'mu.theta', 'sigma.theta')

## run the model
schools.sim <- bugs(data, inits, model.file = model.file1, parameters=parameters,
                    n.chains = 3, n.iter = 1000,
                    ##OpenBUGS.pgm = OpenBUGS.pgm,   ## syntax for R2OpenBUGS::bugs
                    program = 'OpenBUGS',
                    bugs.directory = bugs.directory,
                    WINE = WINE, WINEPATH = WINEPATH, useWINE = TRUE)

## R will pause
## When model is complete a prompt will reappear
print(schools.sim)


## debugging commands in R
## n: next
## c: continue
## s: step into function call
## f: step out of current function call



library(nimble)
modelName <- 'rats'
rats <- readBUGSmodel(model = modelName, dir = getBUGSexampleDir(modelName), returnModelComponentsOnly = TRUE)

code <- rats$model
constants <- rats$data[c('N','T','x')]
data <- rats$data[c('Y')]
inits <- rats$inits

Rmodel <- nimbleModel(code, constants, data, inits, dimensions = rats$dims)

spec <- configureMCMC(Rmodel)
spec$getSamplers()




## creating dendogram figure of autoBlocking for JSM presentation
library(nimble)
load('~/GitHub/legacy/autoBlock/data/model_litters.RData')
Rmodel <- nimbleModel(code, constants, data, inits)
Rmcmc <- buildMCMC(Rmodel, autoBlock = TRUE, makePlots = TRUE)




## working on the conjugacy system, to infer dimensions of 'contributions'
library(nimble)
getDistributionsInfo('dnorm')



## let's manually run all NIMBLE testing
## (to find what's crashing, easier than using travis)
library(nimble)
library(testthat)
setwd('~/GitHub/nimble/nimble/packages/nimble/inst/tests/')
source('test_utils.R')

testFiles <- list.files()
print(testFiles)

##### this listing is out of date!
test_package('nimble', 'copy')      ## pass
test_package('nimble', 'dsl_dists') ## pass
test_package('nimble', 'math')      ## pass
test_package('nimble', 'mcmc')      ## errors/warnings copied below
test_package('nimble', 'meta')      ## pass
test_package('nimble', 'models')    ## 
test_package('nimble', 'trunc')     ## 
test_package('nimble', 'user')      ## 



## fully rebuild NIMBLE user manual

setwd('~/GitHub/nimble/nimble-docs/UserManual/')
library(knitr) 
knit2pdf('NimbleUserManual.Rnw') 
system('open NimbleUserManual.pdf')


## investigating the bug in R M-H sampler

library(nimble)
code <- nimbleCode({
    a ~ dnorm(0, 1)
})
constants <- list()
data <- list()
inits <- list(a = 0)
Rmodel <- nimbleModel(code, constants, data, inits)
spec <- configureMCMC(Rmodel, nodes = NULL)
spec$addSampler('a', 'RW')
Rmcmc <- buildMCMC(spec)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
niter <- 1000

set.seed(0); Rmcmc$run(niter)
set.seed(0); Cmcmc$run(niter)

Rsamples <- as.matrix(Rmcmc$mvSamples)
Csamples <- as.matrix(Cmcmc$mvSamples)

Rsamples
Csamples

d <- Rsamples - Csamples
all(d == 0)



Rmodel$nodes[['a']]$calculateDiff

Rmodel$logProb_a
Cmodel$logProb_a

calculateDiff(Rmodel, 'a')
calculateDiff(Cmodel, 'a')


nimble:::rCalcDiffNodes(Rmodel, 'a')
nimble:::rCalcDiffNodes(Cmodel, 'a')


## testing nimbleOptions() system
## and also the new nimble option: MCMCcontrolDefaultList

library(nimble)

getNimbleOption('MCMCcontrolDefaultList')

environment(sampler_conjugate_dnorm)$nfRefClassDef

getNimbleOption('verifyConjugatePosteriors')
nimble:::setNimbleOption('verifyConjugatePosteriors', TRUE)
getNimbleOption('verifyConjugatePosteriors')
environment(sampler_conjugate_dnorm)$nfRefClassDef
buildConjugateSamplerFunctions()
environment(sampler_conjugate_dnorm)$nfRefClassDef


code <- nimbleCode({
    a ~ dnorm(0, 1)
    b ~ dgamma(a, 2)
    c ~ dbin(b, 10)
    d ~ dnorm(c, 1)
})
constants <- list()
data <- list()
inits <- list(a = 1, b = .5, c = 4, d = 1)

Rmodel <- nimbleModel(code, constants, data, inits)

spec <- configureMCMC(Rmodel)
spec$getSamplers()

lst <- getNimbleOption('MCMCcontrolDefaultList')
lst
lst$sliceWidth <- 99
lst$scale <- 9999
lst
nimbleOptions(MCMCcontrolDefaultList = lst)
getNimbleOption('MCMCcontrolDefaultList')

spec <- configureMCMC(Rmodel)
spec$getSamplers()


## complete example of using RW_llFunction sampler,
## for the User Manual

library(nimble)

code <- nimbleCode({
    p ~ dunif(0, 1)
    y ~ dbin(p, n)
})

Rmodel <- nimbleModel(code, data = list(y = 3),
                      inits = list(p = 0.5, n = 10))

llFun <- nimbleFunction(
    setup = function(model) { },
    run = function() {
        y <- model$y
        p <- model$p
        n <- model$n
        ll <- lfactorial(n) - lfactorial(y) - lfactorial(n-y)
                 + y*log(p) + (n-y)*log(1-p)
        returnType(double())
        return(ll)
    }
)

RllFun <- llFun(Rmodel)

mcmcspec <- configureMCMC(Rmodel, nodes = NULL)

mcmcspec$addSampler(target = 'p', type = 'RW_llFunction',
    control = list(llFunction = RllFun, includesTarget = FALSE))

Rmcmc <- buildMCMC(spec)





## testing new nimble options system

library(nimble)
nimbleOptions()
ls(nimble:::.nimbleOptions)

code <- nimbleCode({
    a ~ dnorm(0, 1)
    b ~ dgamma(a, 1)
    c ~ dbin(b, 10)
    d ~ dnorm(c, 1)
})
constants <- list()
data <- list()
inits <- list(a = 1, b=.5, c=1, d=0)

Rmodel <- nimbleModel(code, constants, data, inits)

spec <- configureMCMC(Rmodel)
spec$getSamplers()
Rmcmc <- buildMCMC(spec)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Cmcmc$run(10000)
samples <- as.matrix(Cmcmc$mvSamples)
apply(samples, 2, mean)



## working on dipper models

rm(list=ls())
library(nimble)
trunc <- TRUE
load('~/GitHub/userDistMCMC/dipperData.RData')
## optionally truncate data:
last <- apply(y, 1, function(hist) max(which(hist==1)))
yDHMM <- 2 - y
onesMatrix = matrix(1,nind,k)  ## 'ones' matrix for survival
onesVector = rep(1,nind)       ## 'ones' vector for chi
if(trunc) { ind <- 1:3;   nind<-length(ind);   first<-first[ind];   last<-last[ind];   y<-y[ind,,drop=FALSE];   yDHMM<-yDHMM[ind,,drop=FALSE];   x_init<-x_init[ind,,drop=FALSE];   onesMatrix<-onesMatrix[ind,,drop=FALSE];   onesVector<-onesVector[ind] }
code <- quote({
    phi ~ dunif(0, 1)
    p ~ dunif(0, 1)
    for(i in 1:nind) {
        for(t in (first[i] + 1):last[i]) {
            y[i,t] ~ dbin(p, 1)
            ##onesMatrix[i,t] ~ dbin(phi, 1)
        }
        ##onesVector[i] ~ dbin(chi[last[i]], 1)
    }
    ##chi[k] <- 1
    ##for(t in 1:(k-1)) {
    ##chi[k-t] <- (1-phi) + phi * (1-p) * chi[k-t+1]
    ##}
})
constants <- list(k=k, nind=nind, first=first, last=last)
data      <- list(y=y)##, onesMatrix=onesMatrix, onesVector=onesVector)
inits     <- list(phi=0.6, p=0.9)

md <- nimbleModel(code, constants, data, inits, returnDef = TRUE)


Rmodel <- nimbleModel(code, constants, data, inits)



## reproducible example of problem in new igraph package
## emailed this example to Gabor


remove.packages('igraph')
install.packages('igraph')

library(nimble)
code <- nimbleCode({
    a ~ dnorm(0, 1)
    b ~ dnorm(a, 1)
})
constants <- list()
data <- list()
inits <- list(a=0, b=0)

Rmodel <- nimbleModel(code, constants, data, inits)


library(igraph)
graph <- graph.empty()
graph <- add.vertices(graph, 2, name = c('a', 'b'))
graph <- add.edges(graph, c(1, 2))
toposortReturn <- topological.sort(graph, mode = 'out')
print(toposortReturn)
class(toposortReturn)



oldGraphID_2_newGraphID <- sort(newGraphID_2_oldGraphID, index = TRUE)$ix
graph <<- permute.vertices(graph, oldGraphID_2_newGraphID)  # re-label vertices in the graph



## writing stan model file for RSBS
##source('http://mc-stan.org/rstan/install.R', echo = TRUE, max.deparse.length = 2000)
##install_rstan()

load('~/GitHub/legacy/automated-blocking-examples/data/model_redblue.RData')
constantsAndData <- c(constants, data)

library(rstan)

stan_mod <- stan_model(file = '~/temp/redblue.stan')

stan_out <- sampling(stan_mod, data=constantsAndData, chains=1, iter=10000, thin=1, init=list(inits))

tempArray <- extract(stan_out, permuted = FALSE, inc_warmup = TRUE)[, 1, ]
dimnames(tempArray)[[2]] <- gsub('_', '.', dimnames(tempArray)[[2]])
if(!all(monitorNodesBUGS %in% dimnames(tempArray)[[2]])) {
    missingNames <- setdiff(monitorNodesBUGS, dimnames(tempArray)[[2]])
    warning(paste0('Stan output is missing values for: ', paste0(missingNames,collapse=', ')))
}
samplesArray <- array(0, dim = c(nkeep, length(monitorNodesBUGS)))
dimnames(samplesArray)[[2]] <- monitorNodesBUGS
monitorsWeHave <- intersect(monitorNodesBUGS, dimnames(tempArray)[[2]])
samplesArray[, monitorsWeHave] <- tempArray[(burnin+1):floor(niter/thin), monitorsWeHave, drop=FALSE]
addToOutput('stan', samplesArray, timeResult)

## playing around to build, install, run, and test gpmanagement package

q('no')
library(devtools)
setwd('~/GitHub/forks/gpmanagement')
remove.packages('gpmanagement')
install_github('danielturek/gpmanagement')
##document('gpmanagement')
system('R CMD build gpmanagement')
system('R CMD INSTALL .')

library(nimble)
library(testthat)
library(gpmanagement)

check('.')  ## this gives me the errors 

gp_setup

test_package('gpmanagement')

source('~/GitHub/forks/gpmanagement/tests/testthat/test_gp.R')

##setwd('~/GitHub/forks/gpmanagement/tests/')
##test_check('gpmanagement')

nf <- gpmanagement:::nfd(3)
nf$run()

set.seed(0)
x <- 1:100
y <- sin(x/5) + rnorm(100, 0.1)
ind <- sort(sample(1:100, 40))
xObs <- x[ind]                ## input
yObs <- y[ind]                ## input
xPred <- c(x, 101:120)        ## input

fit <- gp_setup(xObs, yObs, xPred)


## testing of 3-D array arguments, and their automatic demotion to smaller dim arrays

library(nimble)

nf1 <- nimbleFunction(
    run = function(a = double(2), b = double(3)) {
        print('a11: ', a[1,1])
        print('b111: ', b[1,1,1])
        print(dim(a)[1])
        print(dim(a)[2])
        print(dim(b)[1])
        print(dim(b)[2])
        print(dim(b)[3])
    }
)

nf <- nimbleFunction(
    run = function() {
        a <- nimArray(0, 3, 3)
        declare(b, double(3, c(3,3,3)))
        for(i in 1:3)
            for(j in 1:3) {
                a[i,j] <- 10*i + j
                for(k in 1:3)
                    b[i,j,k] <- 100*i+10*j+k
            }
        print(a)
        ##print(b)
        nf1(a[2,1:2], b[,,])
    }
)

nf()
Cnf <- compileNimble(nf)
Cnf()


## testing distribution dDHMM on a multistate example

CdDHMM <- compileNimble(dDHMM)
CrDHMM <- compileNimble(rDHMM)

checkD <- function(ind) {
    y <- yDHMM
    print(y[ind,])
    x <- y[ind, first[ind]:k]
    length <- length(x)
    phi <- 0.6
    p <- 0.9
    T <- array(c(phi, 1-phi, 0, 1), c(2,2))
    Z <- array(c(p,   1-p,   0, 1), c(2,2))
    pi <- c(1, 0)
    condition <- c(1, 0)
    print(dDHMM(x, length, pi, Z, T, condition, 0))
    print(CdDHMM(x, length, pi, Z, T, condition, 0))
    print(exp(dDHMM(x, length, pi, Z, T, condition, 1)))
    print(exp(CdDHMM(x, length, pi, Z, T, condition, 1)))
}
checkD(225)   ## 0.2484
checkD(202)   ## 0.017496
checkD(167)   ## 0.1246882

checkR <- function(length) {
    phi <- 0.6
    p <- 0.9
    T <- array(c(phi, 1-phi, 0, 1), c(2,2))
    Z <- array(c(p,   1-p,   0, 1), c(2,2))
    pi <- c(1, 0)
    condition <- c(1, 0)
    print(rDHMM(1, length, pi, Z, T, condition))
    print(CrDHMM(1, length, pi, Z, T, condition))
}
checkR(1)
checkR(2)
checkR(3)
checkR(10)


Z <- array(c(.9,.05,.05,.1,.8,.1,0,1/2,1/2), c(3,3))
Z
T <- array(c(0,0,1,0,0,1,0,0,1), c(3,3))
T
pi <- c(1, 1, 1)
condition <- c(1, 1, 1)

y <- c(1,1)
length <- length(y)
dDHMM( y, length, pi, Z, T, condition, 0)
CdDHMM(y, length, pi, Z, T, condition, 0)


## testing new DSL functions:
## nimVector()
## nimArray()

library(nimble)

Rnf <- nimbleFunction(
    run = function(val = double(), len = double(), r = double(), c = double()) {
        onesVec <- nimVector(1, 5)
        print('onesVec:');   print(onesVec)
        zerosVec <- nimVector(0, 5)
        print('zerosVec:');   print(zerosVec)
        otherVec <- nimVector(val, len)
        print('otherVec:');   print(otherVec)
        sumVec <- zerosVec + onesVec + otherVec[1:5]
        print('sumVec:');   print(sumVec)
        arr <- nimArray(.4, r, c)
        print('arr:');   print(arr)
        print('sum of otherVec: ', sum(otherVec))
        print('sum of arrL: ', sum(arr))
    }
)

Cnf <- compileNimble(Rnf)

Rnf(10, 7, 3, 4)
Cnf(10, 7, 3, 4)

## trying nimVector() and nimArray() in BUGS code ?!?
## it works!!
code <- nimbleCode({
    a[1:4] <- nimVector(2, 4)
    phi ~ dnorm(0, 1)
    b[2:10] <- nimVector(10, 9)
    b[1] <- phi
})

Rmodel <- nimbleModel(code)
simulate(Rmodel, 'phi')
calculate(Rmodel)

Cmodel <- compileNimble(Rmodel)

model <- Rmodel
model <- Cmodel

model$a
model$phi
model$b

simulate(model)
calculate(model)

## re-creating NA's data example error
## 11 June 2015:
## Perry said he might get to this by July
## running tests/test-copy.R also gives this same error

library(nimble)

code <- nimbleCode({
    y[2:4] ~ dmnorm(mu3[1:3], cov = cov3[1:3, 1:3])
})
constants <- list(mu3=rep(0,3), cov3=diag(3))
data <- list(y = c(NA,0,0,0))

Rmodel <- nimbleModel(code, constants, data)

calculate(Rmodel)
simulate(Rmodel)

Cmodel <- compileNimble(Rmodel)
calculate(Rmodel)
simulate(Rmodel)


## same error from test-copy.R
library(testthat)
library(nimble)
test_package('nimble', 'copy')




## Zhenglei's example, trying MCMC suite on it

library(lme4)
library(ggplot2)
library(scales)
library(nimble)
library(coda)

I <- 25
J <- 5
beta <- c(-0.2, 0.1)
theta <- 1
sigma <-0.5

simfun <- function(substance=letters[1:10],experiment=factor(1:5),log10HQ=seq(-1,3,by=0.2),theta=1,sigma=0.1,beta=c(-0.2, 0.5)){
    expdat <- expand.grid(experiment = experiment,log10HQ=log10HQ,substance=substance)
    expdat$HQ <- 10^expdat$log10HQ
    expdat$obs <- factor(seq(nrow(expdat)))
    n <- nrow(expdat)
    nsub <- length(substance)
    re <- rnorm(nsub,0,theta)
    names(re) <- substance
    fixed <- beta[1]+beta[2]*expdat$log10HQ
    expdat$y <- re[expdat$substance]+rnorm(n,0,sigma)+fixed
    expdat$effect <- gtools::inv.logit(expdat$y)
    return(expdat)
}

expdat <- simfun(substance=factor(paste0('S',1:I)),experiment = factor(1:J),log10HQ=seq(-1,3,by=0.2),theta=0.5,beta=c(-0.2,0.5),sigma=0.3)

code <- nimbleCode({
    alpha ~ dnorm(0, 0.001)
    beta ~ dnorm( 0, 0.001 )
    tau.y ~ dgamma(0.001, 0.001)  
    tau.re ~ dgamma(0.001, 0.001)
    for(j in 1:M){
        ranef.v[j]~dnorm(0, tau.re)
    }
    for (i in 1:N){
        y[i] ~ dnorm ( mu.y[i], tau.y )
        mu.y[i] <- alpha+beta*x[i] + ranef.v[unit[i]]
    }
})

constants <- list(N = nrow(expdat),M=nlevels(expdat$substance),x = expdat$log10HQ,unit=as.numeric(expdat$substance))
##data <- data.frame(y = expdat$y)
data <- list(y = expdat$y)
inits <- list(alpha = 0, beta = 0.5, tau.y=1,tau.re=1)

out <- MCMCsuite(
    code, constants, data, inits,
    MCMCs = c('nimble', 'jags', 'noConj'),
    niter = 20000,
    makePlot = TRUE,
    savePlot = TRUE,
    summaryStats = c('mean', 'median', 'sd', 'CI95_low', 'CI95_upp', 'effectiveSize'),
    MCMCdefs = list(noConj = quote({ configureMCMC(Rmodel, useConjugacy=FALSE) }))
)

out$summary[, 'effectiveSize', ] / out$timing[c('nimble', 'jags', 'noConj')]


## helping Nick get started with using node$get_mean() function for particle filter

library(nimble)

code <- nimbleCode({
    x ~ dnorm(1, 1)
    y ~ dnorm(x, 1)
})

Rmodel <- nimbleModel(code)
simulate(Rmodel)
Rmodel$x
Rmodel$y
Rmodel$nodeFunctions[['x']]$get_mean()
Rmodel$nodeFunctions[['y']]$get_mean()

## only necessary for the 'higherLevel' example
virtual_NF <- nimbleFunctionVirtual(
    ## don't need to supply prototpye for run(),
    ## since the default is no args, and returnType void
    methods = list(  ## definitely need prototype for return_mean() function
        return_mean = function() {
            returnType(double())
        }
    )
)

nfDef <- nimbleFunction(
    contains = virtual_NF,  ## only necessary for the 'higherLevel' example
    setup = function(model, node) {
        nfList <- nimbleFunctionList(node_stoch_dnorm)
        nfList[[1]] <- model$nodeFunctions[[node]]
    },
    run = function() {
        print('node value is: ', model[[node]])
        print('(assumed to be dnorm) node mean is: ', nfList[[1]]$get_mean())
    },
    methods = list(                        ## only necessary for the 'higherLevel' example
        return_mean = function() {         ##
            returnType(double())           ##
            return(nfList[[1]]$get_mean()) ##
        }                                  ##
    )                                      ##
)

nfX <- nfDef(Rmodel, 'x')
nfY <- nfDef(Rmodel, 'y')
nfX$run()
nfY$run()

compiledList <- compileNimble(list(Rmodel, nfX, nfY))
Cmodel <- compiledList[[1]]
CnfX <- compiledList[[2]]
CnfY <- compiledList[[3]]

CnfX$run()
CnfY$run()

CnfX$return_mean()
CnfY$return_mean()

## this is the 'higherLevel' example
nfDefHigherLevel <- nimbleFunction(
    setup = function(model, nodes) {
        nfList <- nimbleFunctionList(virtual_NF)
        for(i in seq_along(nodes))
            nfList[[i]] <- nfDef(model, nodes[i])
    },
    run = function() {
        declare(meansVector, double(1, length(nfList)))  ## declares a vector of doubles, with length equal to the length of nfList
        for(i in seq_along(nfList)) {
            nfList[[i]]$run()
            meansVector[i] <- nfList[[i]]$return_mean()
        }
        print('vector of the means is: ', meansVector)
    }
)

nfHigherLevel <- nfDefHigherLevel(Rmodel, c('x', 'y'))

nfHigherLevel$run()

CnfHigherLevel <- compileNimble(nfHigherLevel, project = Rmodel)

CnfHigherLevel$run()


## fixing values(model, nodes) <- value reference class '<<-' warning


library(nimble)

nfdef <- nimbleFunction(
    setup = function(model) {},
    run = function(val = double(1)) {
        values(model, 'xxx') <<- val
    }
)

code <- nimbleCode({
    xxx ~ dnorm(0, 1)
})

Rmodel <- nimbleModel(code)

nf <- nfdef(Rmodel)

Cmodel <- compileNimble(Rmodel)
Cnf <- compileNimble(nf, project = Rmodel)

Rmodel$xxx
nf$run(3)
Rmodel$xxx

Cmodel$xxx
Cnf$run(4)
Cmodel$xxx


## testing / trying to find some error in conjugacy system(??)  May 2015

library(nimble)
source('~/GitHub/nimble/packages/nimble/inst/tests/test_utils.R')

## run entire MCMC testing system
source('~/GitHub/nimble/packages/nimble/inst/tests/test-mcmc.R')

set.seed(0)
mu0 = 1:3
Q0 = matrix(c(1, .2, .8, .2, 2, 1, .8, 1, 2), nrow = 3)
Q = solve(matrix(c(3, 1.7, .9, 1.7, 2, .6, .9, .6, 1), nrow = 3))
a = c(-2, .5, 1)
B = matrix(rnorm(9), 3)
code <- nimbleCode({
  mu[1:3] ~ dmnorm(mu0[1:3], Q0[1:3, 1:3])
  y_mean[1:3] <- asCol(a[1:3]) + B[1:3, 1:3] %*% asCol(mu[1:3])
  y[1:3] ~ dmnorm(y_mean[1:3], Q[1:3, 1:3])
})
mu <- mu0 + chol(solve(Q0)) %*% rnorm(3)
# make sure y is a vec not a 1-col matrix or get a dimensionality error
y <- c(a + B%*%mu + chol(solve(Q)) %*% rnorm(3))
data = list(mu0 = mu0, Q0 = Q0, Q = Q, a = a, B = B, y = y)
muQtrue = t(B) %*% Q%*%B + Q0
muMeanTrue = c(solve(muQtrue, crossprod(B, Q%*%(y-a)) + Q0%*%mu0))

constants <- list(mu0 = mu0, Q0 = Q0, Q = Q, B = B, a=a)
data <- list(y=y)
m <- nimbleModel(code, constants = constants, data=data)

m <- nimbleModel(code, constants = data)

spec <- configureMCMC(m)
Rmcmc <- buildMCMC(spec)
spec$getSamplers()

options(error = recover)

test_mcmc(model = code, data = data, seed = 0, numItsC = 10000,
          results = list(mean = list(mu = muMeanTrue),
                           cov = list(mu = solve(muQtrue))),
          resultsTolerance = list(mean = list(mu = rep(.02,3)),
            cov = list(mu = matrix(.01, 3, 3))))


## testing distribution of particle filter likelihood estimates

ll <- log(sapply(muvec, function(mu) mean(dnorm(y, rnorm(n, mu, sigx), sigy))))

plot(muvec, varll)
varest <- 1/2*sigx^4/sigy^4 + muvec^2*sigx^2/sigy^4
lines(muvec, varest, col='red')



## library(nimble)
## Rpf <- buildPF(Rmodel, 'x')
## Cmodel <- compileNimble(Rmodel)
## Cpf <- compileNimble(Rpf, project = Rmodel)
## Rpf$run(10)
## Cpf$run(10)
## ll <- numeric()
## for(i in seq_along(muvec)) {
##     mu <- muvec[i]
##     Cmodel$mu <- mu
##     ll[i] <- Cpf$run(10000)
## }
## plot(muvec, ll)


library(nimble)
source('~/GitHub/pfLL/pfLL.R')
m <- 1000
rep <- 1000
n <- m*rep
muvec <- seq(-5, 5, by=0.2)
sigx <- 1
sigy <- 1
taux <- 1/sigx^2; tauy <- 1/sigy^2
y <- 0   ## observation
constants <- list()
inits <- list(mu = y, sigx = sigx, sigy = sigy)
data <- list(y = y)
code <- quote({
    mu ~ dnorm(0, 0.00001)  ## these should not be necessary
    sigx ~ dunif(0, 1000)  ## these should not be necessary
    sigy ~ dunif(0, 1000)  ## these should not be necessary
    x ~ dnorm(mu, sigx)
    y ~ dnorm(x, sigy)
})
Rmodel <- nimbleModel(code, constants, data, inits)

out <- pfLL(Rmodel, 'x', param = data.frame(mu=muvec), m=m, rep=rep, makePlot=FALSE)

ll <- apply(out$ll, 1, mean)
##ll <- log(sapply(muvec, function(mu) mean(dnorm(y, rnorm(n, mu, sigx), sigy))))
plot(muvec, ll)
ll.pred <- -1/2 * muvec^2 / (sigx^2+sigy^2) + log(1/sqrt(2*pi*(sigx^2+sigy^2)))
lines(muvec, ll.pred, col='red')


varll <- apply(out$ll, 1, var)
##varll <- varll - min(varll)  ### TEMPORARY
plot(muvec, varll)
var.pred <- exp(-8.8 + 1.2592*abs(muvec))
lines(muvec, var.pred, col='red')

1



n <- 100000
sigx <- 1
muvec <- seq(-3, 3, by=0.1)
nmu <- length(muvec)
X <- array(NA, c(n, nmu))
for(i in 1:nmu) X[,i] <- rnorm(n, muvec[i], sigx)

Y <- X^2
var <- apply(Y, 2, var)
pred <- muvec^2 + sigx^2
plot(muvec, var, type='p')
lines(muvec, pred, col='red')



## figuring out 3D plotting in R
rm(list=ls())
df <- expand.grid(list(x=1:10, y=1:10))
rep <- 1
ll <- array(NA, c(dim(df)[1], rep))
for(i in 1:dim(df)[1]) {
    for(j in 1:rep) {
        x <- df[i, 'x']
        y <- df[i, 'y']
        dist <- sqrt((x-4)^2 + (y-7)^2)
        ll[i,j] <- rnorm(1, 100-dist^2, sd = dist/10)
    }
}

library(plot3D)
?scatter3D

x <- rep(df$x, rep)
y <- rep(df$y, rep)
z <- as.numeric(ll)

scatter3D(x, y, z)


## minimally reproducible example of getDependencies error for Perry
library(nimble)
code <- nimbleCode({
    x[1] ~ dnorm(0, 1)
    x[2] ~ dnorm(0, 1)
    y[1] ~ dnorm(x[1], 1)
    y[2] ~ dnorm(x[2], 1)
})

Rmodel <- nimbleModel(code)

Rmodel$getDependencies('x[1]')   ## should be x[1], y[1]
## [1] 'x[1]' 'y[1]' 'y[2]'

Rmodel$getDependencies('x[2]')   ## should be x[2], y[2]
## [1] 'x[2]'



## trying to figure out conjugate sampling in 'equiv' BUGS example model
library(nimble)
code <- nimbleCode({
    tau[1] ~ dgamma(0.001, 0.001)
    tau[2] ~ dgamma(0.001, 0.001)
    pi ~ dnorm(0, 1.0E-06)
    phi ~ dnorm(0, 1.0E-06)
    mu ~ dnorm(0, 1.0E-06)
    for(i in 1:N) { 
        d[i] ~ dnorm(0,tau[2])  # Subject random effect
        for (k in 1:2) {
            Treat[i,k] <- group[i]*(k-1.5) + 1.5  # treatment given
            m[i,k] <- mu + pow(-1, Treat[i,k]-1)* phi /2 + pow(-1, k-1)* pi /2 + d[i]
            Y[i,k] ~ dnorm(m[i,k], tau[1])
        }
    }
    ##theta <- exp(phi)
    ##equivalence <- step(theta - 0.8) - step(theta - 1.2)
})
Y <- structure(c(1.4, 1.64, 1.44, 1.36, 1.65, 1.08, 1.09, 1.25, 1.25, 1.3, 1.65, 1.57, 1.58, 1.68, 1.69, 1.31, 1.43, 1.44, 1.39, 1.52), .Dim = c(10, 2))
group <- c(1, 1, -1, -1, -1, 1, 1, 1, -1, -1)
N <- 10
##T <- 2
constants <- list(N=N, group=group)
data <- list(Y=Y)
inits <- list(mu=0, phi=0, pi=0, tau=c(1,1))
Rmodel <- nimbleModel(code, constants, data, inits)
spec <- configureMCMC(Rmodel)
spec$getSamplers()
Rmcmc <- buildMCMC(spec)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

niter <- 10

set.seed(0); Rmcmc$run(niter)
Rsamples <- as.matrix(Rmcmc$mvSamples)
Rsamples

set.seed(0); Cmcmc$run(niter)
Csamples <- as.matrix(Cmcmc$mvSamples)
Csamples <- Csamples[, dimnames(Rsamples)[[2]]]
Csamples

Rsamples-Csamples


## test equiv model
library(nimble)
source('~/GitHub/nimble/packages/nimble/inst/tests/test_utils.R')
test_mcmc('equiv', numItsC = 1000, resampleData = TRUE)


1

## trying Kalman Filter and calculating likelihood for a linear SSM

# Say state process is
# X(t+1) = r * X(t) + b + nu(t)
# nu(t) ~ dnorm(0, sig2nu)
# 
# Obs process is
# Y(t) = X(t) + eps(t)
# eps(t) ~ dnorm(0, sig2eps)
# 
# Time series goes from t = 1:T
# 
# This is the simplest possible model for us.
# 
# Direct mvn approach:
#   1. Stationary variance of states:
#   V[x] = r^2 V[x] + sig2nu
# V[x] = sig2nu / (1-r^2)
# 
# 2. Deriving entries in covariance matrix of Y (from 1:T)
# Cov[X(t), X(t+1)] = r V[x]
# Cov[X(t), X(t+d)] = r^d V[x]
# Therefore Cov[Y(t), Y(t+d)] = r^d V[x]
# 
# And on the diagonal:
#   V[y] = V[x] + sig2eps
# 
# And mu(t) = E[Y] = E[X] = b/(1-r) for all t
# 
# So with that mean vector and covariance matrix in hand 
# you can directly calculate mvnorm(Y, mu, Cov(Y))
# 
# Kalman filter approach and direct LL are shown in code below.  
# Sorry  I didn't think of this.  I had this from the early 
# explorations that led to Knape, Besbeas and de Valpine in 
# Ecology last year (or was it this year?).  
# 
# In this setup, you want numSites = 1 (we were considering time-series 
# with data from multiple sites).  The LL functions use the meanData 
# (mean over all the sites, which we were comparing to a more complete 
# site-specific model), which for you will be the same as the fullData.



simData <- function(mu=10, a=0.8, sigPN=0.2, sigOE=0.4, t=20) {
    b <- mu*(1-a)
    X <- mu
    for(i in 1:50) X <- a*X + b + rnorm(1,0,sigPN) ## get to stationary
    states <- numeric(t)
    states[1] <- X
    for(i in 2:t) states[i] <- a*states[i-1]+b + rnorm(1,0,sigPN)
    obsStates <- rnorm(t, states, sigOE)
    list(x = states, y = obsStates, mu=mu, a=a, b=b, sigOE=sigOE, sigPN=sigPN, t=t)
}

mvn_ll <- function(d) {
    require(mvtnorm)
    mu <- d$mu
    a <- d$a
    b <- mu*(1-a)
    sigPN <- d$sigPN
    sigPN2 <- sigPN*sigPN
    sigOE <- d$sigOE
    sigOE2 <- sigOE*sigOE
    y <- d$y
    t <- d$t
    mu_vec <- rep( b/(1-a), t)
    cov_mat <- matrix(nrow = t, ncol = t)
    var_x <- sigPN2 / (1-a^2)
    var_y <- var_x + sigOE2
    diag(cov_mat) <- var_y
    for(i in 2:t) {
        for(j in 1:(i-1)) {
            cov_mat[i,j] <- (a^(i-j))*var_x
            cov_mat[j,i] <- cov_mat[i,j]
        }
    }
    dmvnorm(y, mu_vec, cov_mat, log = TRUE)
}

KF_ll <- function(d) {
    mu <- d$mu
    a <- d$a
    b <- mu*(1-a)
    sigPN <- d$sigPN
    sigPN2 <- sigPN*sigPN
    sigOE <- d$sigOE
    sigOE2 <- sigOE*sigOE
    y <- d$y
    t <- d$t
    mu_x <- b/(1-a)
    var_x <- sigPN2 / (1-a^2)
    cov_xy <- var_x
    var_y <- var_x + sigOE2
    ll <- dnorm(mu_x, y[1], sqrt(var_y), log = TRUE)
    mu_x <- mu_x + (cov_xy / var_y) * (y[1] - mu_x)
    var_x <- var_x - cov_xy*cov_xy / var_y
    for(i in 2:t) {
        mu_x <- a*mu_x + b
        var_x <- a^2 * var_x + sigPN2
        if(!is.na(y[i])) {
            cov_xy <- var_x
            var_y <- var_x + sigOE2
            ll <- ll + dnorm(mu_x, y[i], sqrt(var_y), log = TRUE)
            mu_x <- mu_x + (cov_xy / var_y) * (y[i] - mu_x)
            var_x <- var_x - cov_xy*cov_xy / var_y
        }
    }
    ll
}

code <- quote({
    x[1] ~ dnorm(mu, sd = sqrt((sigPN^2) / (1 - a^2)))
    y[1] ~ dnorm(x[1], sd = sigOE)
    for(i in 2:t){
        x[i] ~ dnorm(x[i-1] * a + b, sd = sigPN)
        y[i] ~ dnorm(x[i], sd = sigOE)
    }
})

d <- simData()

library(nimble)
constants <- list(a=d$a, b=d$b, t=d$t, sigOE=d$sigOE, sigPN=d$sigPN, mu=d$mu)
data <- list(y = d$y)
inits <- list(x = d$x)
Rmodel <- nimbleModel(code = code, constants = constants, data=data, inits=inits)
calculate(Rmodel)
pf <- buildPF(Rmodel, 'x')
Cmodel <- compileNimble(Rmodel)
Cpf <- compileNimble(pf, project = Rmodel)

Cpf$run(100000)

mvn_ll(d)

KF_ll(d)



## exploring how truncation works in our system

library(nimble)
code <- nimbleCode({
    a ~ dnorm(0, 1)
    b ~ T(dnorm(a, 1), 1.5, 1.6)
    c ~ T(dnorm(b, 1), 3, Inf)
    d ~ dnorm(c, 1)
    e ~ dnorm(d, 1)
})
constants <- list()
data <- list()
inits <- list()

md <- nimbleModel(code, constants, data, inits, returnDef = TRUE)
Rmodel <- md$newModel(data=data, inits=inits)

Rmodel$modelDef$printDI()

node <- 'c'
Rmodel$nodes[[node]]$calculate
Rmodel$nodes[[node]]$simulate
Rmodel$nodes[[node]]$getLogProb

lapply(Rmodel$modelDef$declInfo, function(di) di$truncation)
for(node in c('a', 'b', 'c')) print(Rmodel$isTruncated(node))
for(node in c('a', 'b', 'c')) print(Rmodel$getBounds(node))

spec <- configureMCMC(Rmodel)
spec$getSamplers()
Rmcmc <- buildMCMC(spec)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Cmcmc$run(10000)
samples <- as.matrix(Cmcmc$mvSamples)
apply(samples, 2, mean)



## using rankSample

library(nimble)
##Cmodel <- compileNimble(Rmodel = nimbleModel(quote({a ~ dnorm(0, 1)})))
nfDef <- nimbleFunction(
    setup = function() {},
    run = function(wts = double(1), m = integer(), silent = logical()) {
        ##print('in nimbleFunction rankSample')
        declare(samp, integer(1))
        rankSample(wts, m, samp, silent)
        print(samp)
    }
)
nf <- nfDef()
Cnf <- compileNimble(nf)

nfDefOmit <- nimbleFunction(
    setup = function() {},
    run = function(wts = double(1), m = integer()) {
        ##print('in nimbleFunction rankSample')
        declare(samp, integer(1))
        rankSample(wts, m, samp)
        print(samp)
    }
)
nfOmit <- nfDefOmit()
CnfOmit <- compileNimble(nfOmit)
x <- 1:10

silent <- FALSE
silent <- TRUE

set.seed(1);   .Call('rankSample', c(.5,.4,.1), 10L, x, silent)
set.seed(1);            rankSample(c(.5,.4,.1), 10L, x, silent);    print(x)
set.seed(1);                nf$run(c(.5,.4,.1), 10L,    silent)
set.seed(1);               Cnf$run(c(.5,.4,.1), 10L,    silent)

set.seed(1);   .Call('rankSample', c( 0, 0, 0), 10L, x, silent)
set.seed(1);            rankSample(c( 0, 0, 0), 10L, x, silent);    print(x)
set.seed(1);                nf$run(c( 0, 0, 0), 10L,    silent)
set.seed(1);               Cnf$run(c( 0, 0, 0), 10L,    silent)

set.seed(1);   .Call('rankSample', c(.5,.4,-1), 10L, x, silent)
set.seed(1);            rankSample(c(.5,.4,-1), 10L, x, silent);    print(x)
set.seed(1);                nf$run(c(.5,.4,-1), 10L,    silent)
set.seed(1);               Cnf$run(c(.5,.4,-1), 10L,    silent)

set.seed(1);   .Call('rankSample', c(.5,.4,-1), 10L, x)
set.seed(1);            rankSample(c(.5,.4,-1), 10L, x);            print(x)
set.seed(1);            nfOmit$run(c(.5,.4,-1), 10L)
set.seed(1);           CnfOmit$run(c(.5,.4,-1), 10L)



## difference in rank sample
library(nimble)
nfDef <- nimbleFunction(
    setup = function() {},
    run = function() {
        declare(wts, double(1, 1))
        wts[1] <- 1
        m <- 1L
        declare(samp, integer(1, 1))
        rand <- rnorm(1, 0, 1)
        print('random number before rankSample: ', rand)
        rankSample(wts, m, samp)
        rand <- rnorm(1, 0, 1)
        print('random number after rankSample: ', rand)
    }
)

nf <- nfDef()
Cnf <- compileNimble(nf)
set.seed(0)
nf$run()
set.seed(0)
Cnf$run()



## testing buildPF() particle filter (pf) algorithm
library(nimble)
code <- quote({
    mu ~ dnorm(0, sd = 1000)
    b ~ dnorm(0, sd = 1000)
    sigPN ~ dunif(0.0001, 1)
    sigOE ~ dunif(0.0001, 1)
    x[1] ~ dnorm(mu, sd = sqrt(sigPN^2 + sigOE^2))
    y[1] ~ dnorm(x[1], sd = sigOE)
    a <- 1-(b/mu)
    for(i in 2:t){
        x[i] ~ dnorm(x[i-1] * a + b, sd = sigPN)
        y[i] ~ dnorm(x[i], sd = sigOE)
    }
})
t <- 5
constants <- list(t = t)
Rmodel <- nimbleModel(code = code, constants = constants)
Rmodel$mu <- 1/(1-.95)
Rmodel$b <- 1
Rmodel$sigPN <- .2
Rmodel$sigOE <- .05
set.seed(0)
calculate(Rmodel, Rmodel$getDependencies(c('mu','b','sigPN','sigOE'), determOnly = TRUE))
simulate(Rmodel, Rmodel$getDependencies(c('x', 'y')))
data <- list(y = Rmodel$y)
inits <- list(mu = Rmodel$mu, b = Rmodel$b, sigPN = Rmodel$sigPN, sigOE = Rmodel$sigOE, x = Rmodel$x)
rm(Rmodel)
Rmodel <- nimbleModel(code, constants, data, inits)
pf <- buildPF(Rmodel, 'x')
Cmodel <- compileNimble(Rmodel)
Cpf <- compileNimble(pf, project = Rmodel)

m <- 100
set.seed(0);    pf$run(m)
set.seed(0);   Cpf$run(m)

1
2
3



## making MCMCsuite work for Ryan using Stan MCMC

allModels <- c('blocker', 'bones', 'dyes', 'line', 'pump', 'rats')
library(nimble)
library(coda)
allModels <- 'blocker'
mcmcs=c('nimble','jags','noConj','stan')
x=list()
for (i in 1:length(allModels)){
    x[[i]]<-readBUGSmodel(model=allModels[i],
                          dir=getBUGSexampleDir(allModels[i]),
                          returnModelComponentsOnly=TRUE)
}
i <- 1
suite_output <- MCMCsuite(
    x[[i]]$model,
    constants = x[[i]]$data,
    inits = x[[i]]$inits,
    MCMCs = mcmcs,
    makePlot=F,
    savePlot=F,
    summaryStats=c('mean','median','sd','CI95_low','CI95_upp','effectiveSize'),
    MCMCdefs = list(noConj = quote({ configureMCMC(Rmodel, useConjugacy=FALSE) })),
    stan_model = paste('~/temp/', allModels[i], '.stan', sep='')
)


## learning about R's memory management and storage of objects
gcinfo(TRUE)
gc()

library(pryr)
?object_size
object_size(1:100)
object_size(mtcars)

x <- 0:50
s <- sapply(x, function(i) object_size(1:i))
plot(x, s, type='s', ylim=c(0,300))

obj <- complex()
object_size(obj)

mem_used
mem_used()

pryr:::node_size()
pryr:::show_bytes

mem_change(NA)

bytes(1)
bytes(1L)
bytes('ab')

address(2)
x <- 1:10
address(x)
refs(x)
y <- x
refs(x)
z <- x
x[1] <- 16L

tracemem(x)
x[1] <- .1

x <- 15
address(x)

system(paste0('~/temp/t ', address(x)))



## testing of mcmcplots package function mcmcplot()

rm(list=ls())
library(nimble)
model <- 'SSMcorrelated'
load(paste0('~/GitHub/autoBlock/data/model_', model, '.RData'))
Rmodel <- nimbleModel(code, constants, data, inits)

out <- MCMCsuite(code, constants, data, inits, monitors=c('a','b','p'), MCMCs=c('nimble','nimble_RW', 'nimble_slice', 'jags'), makePlot=FALSE, niter=10000)

library(coda)
samples <- out$samples
dim(samples)
dimnames(samples)
mcmcList <- list()
for(i in 1:dim(samples)[1])     mcmcList[[i]] <- coda::mcmc(t(samples[i,,]))
names(mcmcList) <- dimnames(samples)[[1]]
codaList <- coda::as.mcmc.list(mcmcList)

library(mcmcplots)
mcmcplot(codaList)
traplot(codaList)




## testing conjugacy in newNimbleModel branch

library(nimble)
code <- nimbleCode({
    for(i in 1:10) {
        a[i] ~ dnorm(0, 1)
        b[i] ~ dgamma(1, 1)
    }
    for(i in 1:5)   a_temp[i] <- 3*a[i]
    for(i in 6:10)  a_temp[i] <- a[i]^2
    for(i in 1:10)  x[i] ~ dnorm(a_temp[i], 1)
    for(i in 1:10) {
        y1[i] ~ dpois(2*b[i])
        y2[i] ~ dgamma(1, rate = b[i])
    }
})
Rmodel <- nimbleModel(code)

spec <- configureMCMC(Rmodel)

spec$getSamplers()
for(i in 1:200) spec$addSampler('RW', 'x[2]', print=F)
spec$addSampler('RW', 'x[2]')

spec$setSamplers(c(1:40, 1:40))
spec$setSamplers(c('a', 'y2'))
spec$removeSamplers('a')
spec$setSamplers()


## NEED TO RE-RUN THIS TEST, ONCE PERRY FIXES THE MODELVALUES COPYING ISSUE
## Dirichlet-multinomial conjugacy
## single multinomial
library(nimble)
set.seed(0)
n <- 100
alpha <- c(10, 30, 15, 60, 1)
K <- length(alpha)
p <- c(.12, .24, .09, .54, .01)
y <- rmulti(1, n, p)
code <- nimbleCode({
    y[1:K] ~ dmulti(p[1:K], n)
    p[1:K] ~ ddirch(alpha[1:K])
    for(i in 1:K) {
        alpha[i] ~ dgamma(.001, .001)
    }
})
constants <- list(n = n, K = K)
data <- list(y = y)
inits <- list(p = rep(1/K, K), alpha = rep(K, K))
Rmodel <- nimbleModel(code, constants=constants, data=data, inits=inits)

mv <- modelValues(Rmodel)
Rmodel$y
mv$y
nimCopy(from=Rmodel, to=mv, nodes = 'y', row = 1)
mv$y


spec <- configureMCMC(Rmodel, monitors = c('alpha', 'p'))
spec$getSamplers()
Rmcmc <- buildMCMC(spec)

set.seed(0)
Rmcmc$run(10)

conjugate posterior density appears to be wrong, off by Inf
conjugate posterior density appears to be wrong, off by Inf
There were 50 or more warnings (use warnings() to see the first 50)

test_mcmc(model = code, data= data, seed = 0, numItsC = 10000,
          inits = inits,
          results = list(mean = list(p = p)),
          resultsTolerance = list(mean = list(p = rep(.06, K))))




## just a good old fashioned test that it's working
library(nimble)
code <- BUGScode({
    x ~ dgamma(1, 1)       # should satisfy 'gamma' conjugacy class
    a  ~ dnorm(0, x)     # should satisfy 'norm' conjugacy class
    a2 ~ dnorm(0, tau = 3*x+0)
    b  ~ dpois(0+5*x)
    b2 ~ dpois(1*x*1)
    c ~ dgamma(1, 7*x*5)
    for(i in 2:3) {
        jTau[i] <- 1
        jNorm[i] ~ dnorm(c * (a+3) - i, var = jTau[i])
        kTauSd[i] <- 2
        kLogNorm[i] ~ dlnorm(0 - a - 6*i, kTauSd[i])
    }
})

constants <- list()
data <- list()
inits <- list()
Rmodel <- nimbleModel(code, constants, data, inits)
Rmodel$checkConjugacy()

spec <- configureMCMC(Rmodel, control = list(scale=0.01), monitors = c('x', 'c'))
spec$getSamplers()
Rmcmc <- buildMCMC(spec)
set.seed(0)
Rmcmc$run(10)
samples <- as.matrix(Rmcmc$mvSamples)
samples


sampleVals = list(x = c(3.950556165467749, 1.556947815895538, 1.598959152023738, 2.223758981790340, 2.386291653164086, 3.266282048060261, 3.064019155073057, 3.229661999356182, 1.985990552839427, 2.057249437940977),
  c = c( 0.010341199485849559, 0.010341199485849559, 0.003846483017887228, 0.003846483017887228, 0.007257679932131476, 0.009680314740728335, 0.012594777095902964, 0.012594777095902964, 0.018179641351556003, 0.018179641351556003))

test_mcmc(model = code, exactSample = sampleVals, seed = 0, mcmcControl = list(scale=0.01))




## testing if a Ref Class can just return a smaller object (not the whole Ref Class object)
RC <- setRefClass(
    Class = 'RC',
    fields = list(a='ANY'),
    methods = list(
        initialize = function(a) {
            a <<- a
            return(a)
        }
    )
)

rc <- RC(1)
class(rc)
rc$a



## bunch of random stuff

fileToList <- function(file) {
    env <- new.env()
    source(fileName, local = env)
    lst <- list()
    for(name in ls(env))   lst[[name]] <- get(name, env)
    return(lst)
}

library(nimble)

NF <- nimbleFunction(
    setup = function() { },
    run = function(a = double(0)) {
        returnType(double(0))
        return(a)
    }
)

myNF <- NF()
myNF_C <- compileNimble(myNF)

> myNF$run
## function (a) 
##     return(a)

myNF_C$run
## function (a) 
## {
##     ans <- .Call('CALL_nfRefClass60_operator_', a, .basePtr)
##     ans <- ans[[2]]
##     ans
## }


### testing the new autoBlock
rm(list=ls())
library(nimble)
model <- 'litters'
load(paste0('~/GitHub/autoBlock/data/model_', model, '.RData'))
Rmodel <- nimbleModel(code, constants, data, inits)
Rmodel$getNodeNames()

ab <- autoBlock(Rmodel, run = runList)
ab
ab$spec$getSamplers()
Rmcmc <- buildMCMC(ab$spec)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project=Rmodel)


spec <- configureMCMC(Rmodel, autoBlock=TRUE)
spec$getSamplers()
Rmcmc <- buildMCMC(spec)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project=Rmodel)

Rmcmc <- buildMCMC(Rmodel, autoBlock=TRUE)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project=Rmodel)







## March 2015
## timing comparisons for red-state-blue-state
library(nimble)
library(rjags)
rm(list=ls())
load('~/GitHub/autoBlock/data/model_redblue.RData')

## this version of red-state-blue-state code is 'jags friendly'
## reparametrized to use Precision matrix, replaced expit with logit
code <- nimbleCode({
    for (i in 1:2) {
        for (j in 1:2) {
            gamma[i, j] ~ dnorm(0, 1e-04)
        }
    }
    sigmaIntercept ~ dunif(0, 100)
    sigmaSlope ~ dunif(0, 100)
    rho ~ dunif(-1, 1)
    Precision[1, 1] <- 1/(1-rho^2) * 1/(sigmaIntercept^2)
    Precision[2, 2] <- 1/(1-rho^2) * 1/(sigmaSlope^2)
    Precision[1, 2] <- 1/(1-rho^2) * -rho/(sigmaIntercept*sigmaSlope)
    Precision[2, 1] <- 1/(1-rho^2) * -rho/(sigmaIntercept*sigmaSlope)
    for (i in 1:Nstates) {
        stateBetaMeans[i, 1] <- gamma[1, 1] + gamma[1, 2] * stateIncome[i]
        stateBetaMeans[i, 2] <- gamma[2, 1] + gamma[2, 2] * stateIncome[i]
        stateBetas[i, 1:2] ~ dmnorm(stateBetaMeans[i, 1:2], Precision[1:2, 1:2])
    }
    for (i in 1:N) {
        logit(p[i]) <- stateBetas[state[i], 1] + stateBetas[state[i], 2] * income[i]
        y[i] ~ dbern(p[i])
    }
})

niter <- 50000
monitorVars <- c('sigmaIntercept', 'sigmaSlope', 'rho', 'gamma')
constsAndData <- c(constants, data)
modelfile <- file.path(tempdir(), 'model.txt')
writeLines(paste0('model\n', paste0(deparse(code, width.cutoff=500L), collapse='\n')), con=modelfile)


system.time(modelDef <- nimbleModel(code=code, constants=constants, data=data, inits=inits, returnDef=TRUE))[3]/60
26

system.time(Rmodel <- modelDef$newModel(data=data, inits=inits))[3]/60
5

system.time(spec <- configureMCMC(Rmodel))
6 seconds

system.time(Rmcmc <- buildMCMC(spec))[3]
15 seconds

system.time(Cmodel <- compileNimble(Rmodel))[3]/60
11

system.time(Cmcmc <- compileNimble(Rmcmc, project = Rmodel))[3]/60
1

system.time(Cmcmc$run(niter))[3]/60
11 minutes for 50,000 iterations

system.time(jags_mod <- jags.model(file=modelfile, data=constsAndData, inits=inits, n.chains=1, quiet=FALSE))[3]/60
10 seconds

system.time(jags_out <- coda.samples(model=jags_mod, variable.names=monitorVars, n.iter=niter, thin=1))[3]/60
2.001833





## displaying calc and sim functions
code <- nimbleCode({ x ~ dnorm(3, sd = 5) })
model <- nimbleModel(code)
model$nodes$x$simulate
## function()
##     model$x <<- rnorm(1, mean = 3, sd = 5)
model$nodes$x$calculate
## function() {
##     model$logProb_x <<- dnorm(model$x, mean = 3, sd = 5, log = 1)
##     return(invisible(model$logProb_x))
## }
model$nodes$x$getLogProb
## function() 
##     return(model$logProb_x)


### March 2015
### trying to figure out new MCMCspec

library(nimble)


#These first few lines are nothing new
myModelCode <- nimbleCode({
	a[1] ~ dnorm(0,1)
	a[2] ~ dexp(a[1]^2)
	a[3] ~ dnorm(0,1)
})
myModel <- nimbleModel(myModelCode)
mySpec1 <- configureMCMC(myModel)
myMCMC1 <- buildMCMC(mySpec1)
cmod <- compileNimble(myModel)
cmcmc1 <- compileNimble(myMCMC1, project = myModel)


# Rebuilding the way we are supposed to: 
# building a new spec with configureMCMC 
# adding a sampler which has already been compiled 
# not altering the monitors
mySpec2 <- configureMCMC(oldSpec = mySpec1)
mySpec2$addSampler('RW', control = list(targetNode = 'a[3]'))
rmcmc2 <- buildMCMC(mySpec2)
cmcmc2 <- compileNimble(rmcmc2, project = myModel)




#### March 2015 working through a full compileNimble call

library(nimble)

code <- nimbleCode({
    a ~ dnorm(0, 1)
})
inits <- list(a = 1)
Rmodel <- nimbleModel(code, inits = inits)

debug(compileNimble)
Cmodel <- compileNimble(Rmodel)

debug(project$compileModel)
debug(modelCpp$buildAll)
debug(buildNodes)
debug(nimbleProject$compileNimbleFunctionMulti)
debug(compileNimbleFunction)
debug(buildNimbleFunctionCompilationInfo)


#### March 2015 testing of RW_llFunction or Carl B, submitted as test to testing suite

library(nimble)

code <- nimbleCode({
    a ~ dnorm(0, 1)
})
inits <- list(a = 0)
Rmodel <- nimbleModel(code=code, inits=inits)

llFunction <- nimbleFunction(
    setup = function(model) { },
    run = function() {
        ll <- dnorm(1, model$a, 1, log=1)
        returnType(double())
        return(ll)
    }
)
 
myLL <- llFunction(Rmodel)

spec <- configureMCMC(Rmodel, nodes=NULL)
spec$addSampler('RW_llFunction', list(targetNode='a', llFunction=myLL, includesTarget=FALSE))
Rmcmc <- buildMCMC(spec)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Rmcmc$run(10)
unlist(Rmcmc$mvSamples[['a']])
## [1] -0.2849616  0.4851313  0.1241257  0.1241257  0.1241257 -0.3553077 0.4223091  0.4455515  0.4455515  0.4455515

set.seed(0)
Cmcmc$run(10)
unlist(Cmcmc$mvSamples[['a']])
## [1] -0.2849616  0.4851313  0.1241257  0.1241257  0.1241257 -0.3553077 0.4223091  0.4455515  0.4455515  0.4455515

Cmcmc$run(100000)
samples <- as.matrix(Cmcmc$mvSamples)
apply(samples, 2, mean)
## should be within some tolerance of 0.5


###### Feb 2015 example for Cal Poly interview talk

library(nimble)

code <- nimbleCode({
    for( i in 1 : N ) {
        b[i] ~ dnorm(mu, tau)
        r[i] ~ dbin(p[i], n[i])
        logit(p[i]) <- b[i]
    }
    pop.mean <- 1 / (1 + exp(-mu))
    mu ~ dnorm(0, 0.001)
    sigma <- 1 / sqrt(tau)
    tau ~ dgamma(0.001, 0.001)
})

constants <- list(N = 10, n = c(12,20,18,19,10,23,19,15,8,10))
data <- list(r = c(4,7,10,3,6,7,11,9,3,8))
inits <- list(mu = 0, tau = 1)

model <- nimbleModel(code=code, constants=constants, data=data, inits=inits)

spec <- configureMCMC(model)
spec$getSamplers()
spec$addSampler(‘slice’, list(targetNode = ‘mu’))
spec$getSamplers()
mcmc <- buildMCMC(spec)

Cmodel <- compileNimble(model)
Cmcmc <- compileNimble(mcmc, project = model)

Cmcmc$run(10000)

samples <- as.matrix(Cmcmc$mvSamples)

apply(samples, 2, mean)

