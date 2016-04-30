
niter <- 10000
monitorVars <- c('a', 'x')

constsAndData <- c(constants, data)
modelfile <- file.path(tempdir(), 'model.txt')
writeLines(paste0('model\n', paste0(deparse(code, width.cutoff=500L), collapse='\n')), con=modelfile)

library(rjags)
jags_mod <- jags.model(file=modelfile, data=constsAndData, inits=inits, n.chains=1, quiet=FALSE)
jags_out <- coda.samples(model=jags_mod, variable.names=monitorVars, n.iter=niter, thin=1)

dimnames(jags_out[[1]])
means <- apply(jags_out[[1]][,], 2, mean)
means
sds <- apply(jags_out[[1]][,], 2, sd)
sds
