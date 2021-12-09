cd ~/github/nimble/nimbleHMC/

rm -f nimbleHMC_*.tar.gz

R CMD build nimbleHMC

R CMD INSTALL nimbleHMC_*.tar.gz -l ~/Documents/

## use in R scripts:
## if(Sys.info()['nodename'] == 'gandalf') library(nimbleHMC, lib.loc = '~/Documents/') else library(nimbleHMC)

