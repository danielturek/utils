cd ~/github/nimble/nimble/

git pull

cd ~/github/nimble/nimble/packages/

rm -f nimble/inst/CppCode/*.o
rm -f nimble/inst/CppCode/*.so
rm -f nimble/src/*.o
rm -f nimble/src/*.so

R CMD build nimble
R CMD INSTALL nimble -l ~/Documents/

## use in R scripts:
## if(Sys.info()['nodename'] == 'gandalf') library(nimble, lib.loc = '~/Documents/') else library(nimble)


## MAKING NEW MCMC GOLD FILE ON GANDALF (A) and (B) below.
##
## (A) change frontmatter of test-mcmc.R, to have it run ond gandalf,
## and write new MCMC gold (test comparison output) file to 'mcmcTestLog_Correct_NEW.Rout' :

#### ##[this is the top of test-mcmc.R]
#### ##source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))
#### library(nimble, lib.loc = '~/Documents/')
#### source('~/github/nimble/nimble/packages/nimble/tests/testthat/test_utils.R')
####  
#### context("Testing of default MCMC")
####  
#### .......
#### .......
#### .......
####  
#### ## If you do *not* want to write to results files
#### ##    comment out the sink() call below.  And consider setting verbose = FALSE 
#### ## To record a new gold file, nimbleOptions('generateGoldFileForMCMCtesting') should contain the path to the directory where you want to put it
#### ## e.g. nimbleOptions(generateGoldFileForMCMCtesting = getwd())
#### ## Comparison to the gold file won't work until it is installed with the package.
####  
#### nimbleOptions(generateGoldFileForMCMCtesting = '~/github/nimble/nimble/packages/nimble/tests/testthat')
####  
#### goldFileName <- 'mcmcTestLog_Correct_NEW.Rout'
#### tempFileName <- 'mcmcTestLog.Rout'
#### generatingGoldFile <- !is.null(nimbleOptions('generateGoldFileForMCMCtesting'))
#### outputFile <- if(generatingGoldFile) file.path(nimbleOptions('generateGoldFileForMCMCtesting'), goldFileName) else tempFileName

## (B) locally on my machine, edit the 'mcmcTestLog_Correct_NEW.Rout' file,
## to remove the 'Test passed' outputs from test that.
## below code is from scratch.R, to accomplish that:

#### ## update a new MCMC testing gold file
#### ## to remove the "Test passed" messages from test_that testthat package
#### ## to compare gold file
#### path <- '~/github/nimble/nimble/packages/nimble/tests/testthat/'
#### infile <- 'mcmcTestLog_Correct_NEW.Rout'
#### outfile <- 'mcmcTestLog_Correct_NEW2.Rout'
#### t <- readLines(paste0(path, infile))
####  
#### i <- 1
#### tot <- length(t)
#### while(i <= tot) {
####     if(grepl('Test passed', t[i])) {
####         t[i] <- gsub('Test passed .', '', t[i])
####         if(i < tot) {
####             t[i] <- paste0(t[i], t[i+1])
####             t <- t[-(i+1)]
####         } else {
####             if(t[i] == '') t <- t[-i]
####         }
####         tot <- length(t)
####     } else {
####         i <- i+1
####     }
#### }
####  
#### grep('Test passed', t, value = TRUE)
#### writeLines(t, con = paste0(path, outfile))
