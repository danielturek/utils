
niter <- 10000
monitorVars <- c('a', 'x')

constsAndData <- c(constants, data)
modelfile <- file.path(tempdir(), 'model.txt')
writeLines(paste0('model\n', paste0(deparse(code, width.cutoff=500L), collapse='\n')), con=modelfile)

library(rjags)
jags_mod <- jags.model(file=modelfile, data=constsAndData, inits=inits, n.chains=1, quiet=FALSE)
list.samplers(jags_mod)
jags_out <- coda.samples(model=jags_mod, variable.names=monitorVars, n.iter=niter, thin=1)

samples <- as.matrix(jags_out[[1]])
dimnames(samples)
dim(samples)
head(samples)

means <- apply(samples, 2, mean)
means
sds <- apply(samples, 2, sd)
sds

## alternate, using jagsUI::jags
library(jagsUI)
set.seed(0)
out <- jags(
    model.file = modelfile,
    data = constsAndData,
    inits = list(inits),   ## a list, of length = n.chains
    parameters.to.save = params,
    n.iter = niter,
    n.chains = 1,
    n.burnin = 0,
    n.thin = 1)
print(out, dig = 3)
## out$sims.list is a *named list* of samples for each parameter
