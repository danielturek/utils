library(nimble)
library(basicMCMCplots)
library(coda)

code <- nimbleCode({
    a ~ dnorm(0, 1)
})
constants <- list()
data <- list()
inits <- list(a = 0)

Rmodel <- nimbleModel(code, constants, data, inits)
Rmodel$calculate()

conf <- configureMCMC(Rmodel)
conf$printSamplers()
conf$printSamplers(byType = TRUE)
conf$printMonitors()

Rmcmc <- buildMCMC(conf)

compiledList <- compileNimble(list(model=Rmodel, mcmc=Rmcmc))
Cmodel <- compiledList$model; Cmcmc <- compiledList$mcmc
##Cmodel <- compileNimble(Rmodel)
##Cmcmc <- compileNimble(Rmcmc, project = Rmodel)#, showCompilerOutput = TRUE)

set.seed(0)
samples <- runMCMC(Cmcmc, 10000)

colnames(samples)
samplesSummary(samples)
samplesPlot(samples)
apply(samples, 2, effectiveSize)


nfDef <- nimbleFunction(
    setup = function() {},
    run = function() {
        returnType()
    }
)

Rnf <- nfDef()
Cnf <- compileNimble(Rnf)#, showCompilerOutput = TRUE)

Rnf$run()
Cnf$run()


Rnf <- nimbleFunction(
    run = function() {
        returnType()
    }
)

Cnf <- compileNimble(Rnf)#, showCompilerOutput = TRUE)

Rnf()
Cnf()

