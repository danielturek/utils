
library(nimble)
library(coda)

niter <- 10000
burnin <- 20000

out1 <- compareMCMCs(modelInfo1, MCMCs = 'nimble', monitors = params,
                     niter = niter, burnin = burnin, summary = FALSE)

out1[[1]] <- rename_MCMC_comparison_method('nimble', 'MCMC1', out1[[1]])

out2 <- compareMCMCs(modelInfo2, MCMCs = 'nimble', monitors = params,
                     niter = niter, burnin = burnin, summary = FALSE)

out2[[1]] <- rename_MCMC_comparison_method('nimble', 'MCMC2', out2[[1]])

combined <- combine_MCMC_comparison_results(out1[[1]], out2[[1]])

make_MCMC_comparison_pages(combined, dir = '~/temp', modelNames = 'modelName')

browseURL(Sys.glob(file.path('modelName.html')))


## ==============================================
## ==============================================

outList <- list()
niter <- 10000



name <- 'default'
out <- compareMCMCs(modelInfo = modelInfo, MCMCs = 'nimble', niter = niter)[[1]]
out <- rename_MCMC_comparison_method('nimble', name, out)
outList[[name]] <- out

name <- 'MH'
out <- compareMCMCs(modelInfo = modelInfo, MCMCs = 'MH', niter = niter,
                    MCMCdefs = list(MH = quote({
                        conf <- configureMCMC(Rmodel, useConjugacy = FALSE)
                        conf$printSamplers()
                        return(conf)})))[[1]]
outList[[name]] <- out


name <- 'langevin'
out <- compareMCMCs(modelInfo = modelInfo, MCMCs = 'langevin', niter = niter,
                    MCMCdefs = list(langevin = quote({
                        conf <- configureMCMC(Rmodel, nodes = NULL)
                        nodeNames <- Rmodel$getNodeNames(stochOnly = TRUE, includeData = FALSE)
                        for(nn in nodeNames) conf$addSampler(nn, 'langevin2')
                        conf$printSamplers()
                        return(conf)})))[[1]]
outList[[name]] <- out

name <- 'langTop'
out <- compareMCMCs(modelInfo = modelInfo, MCMCs = 'langTop', niter = niter,
                    MCMCdefs = list(langTop = quote({
                        conf <- configureMCMC(Rmodel)
                        nodeNames <- c('tau.c', 'tau.beta', 'tau.alpha', 'beta.c', 'alpha.c')
                        conf$removeSamplers(nodeNames)
                        for(nn in nodeNames) conf$addSampler(nn, 'langevin2')
                        conf$printSamplers()
                        return(conf)})))[[1]]
outList[[name]] <- out

results <- do.call(combine_MCMC_comparison_results, unname(outList))

make_MCMC_comparison_pages(results, pageComponents = list(timing = TRUE, efficiencySummary = FALSE, efficiencySummaryAllParams = TRUE, paceSummaryAllParams = TRUE, efficiencyDetails = TRUE, posteriorSummary = TRUE))

system(paste0('open MCMCresults.html'))




##runComparison <- function(modelInfoFile, reduced, name, MCMCs, niter, MCMCdefs = list(), add = FALSE, saveFile, verbose = TRUE) {
##    if(length(MCMCs) > 1) stop('only one MCMC at a time, please')
##    if(verbose) message(paste0('running ', name, ' on ', modelInfoFile, '...'))
##    modelInfoFileToLoad <- modelInfoFile
##    if(reduced) modelInfoFileToLoad <- paste0(modelInfoFileToLoad, '_reduced')
##    modelInfoFileToLoad <- paste0('data/modelInfo_', modelInfoFileToLoad, '.RData')
##    load(modelInfoFileToLoad)
##    outList <- if(add) dget(saveFile) else list()
##    out <- compareMCMCs(modelInfo = modelInfo, MCMCs = MCMCs, MCMCdefs = MCMCdefs,
##                        monitors = modelInfo$monitors, niter = niter)[[1]]
##    out <- rename_MCMC_comparison_method(MCMCs, name, out)
##    outList[[name]] <- out
##    if(!missing(saveFile)) dput(outList, file = saveFile)
##    if(verbose) message(paste0('finished running ', name, ' on ', modelInfoFile))
##    return(invisible(outList))
##}
## 
##makePages <- function(saveFile, dir, open = TRUE) {
##    outList <- dget(saveFile)
##    results <- do.call(combine_MCMC_comparison_results, unname(outList))
##    pagesDir <- paste0('pages/', dir, '/')
##    make_MCMC_comparison_pages(results, dir = pagesDir, pageComponents = list(timing = TRUE, efficiencySummary = FALSE, efficiencySummaryAllParams = TRUE, paceSummaryAllParams = TRUE, efficiencyDetails = TRUE, posteriorSummary = TRUE))
##    if(open) system(paste0('open ', pagesDir, 'MCMCresults.html'))
##}







