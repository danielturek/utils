library(nimble)

code <- nimbleCode({
    a ~ dnorm(0, 1)
})
constants <- list()
data <- list()
inits <- list(a = 0)

Rmodel <- nimbleModel(code, constants, data, inits)

conf <- configureMCMC(Rmodel)
conf$printSamplers()
Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
samples <- runMCMC(Cmcmc, 10000)
##Cmcmc$run(10000)
##samples <- as.matrix(Cmcmc$mvSamples)

colnames(samples)
apply(samples, 2, mean)

samplesPlot(samples)

library(coda)
apply(samples, 2, effectiveSize)

