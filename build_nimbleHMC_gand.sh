
rm -rf ~/Documents/nimble

cd ~/github/nimble/nimble/packages
git checkout ADoak_without_HMC
git pull
rm -f nimble_*.tar.gz
R CMD build nimble
R CMD INSTALL -l ~/Documents nimble

cd ~/github/nimble/nimbleHMC
git checkout master
git pull
rm -f nimbleHMC_*.tar.gz
R CMD build nimbleHMC
R CMD INSTALL -l ~/Documents nimbleHMC

## use in R scripts:
## if(Sys.info()['nodename'] == 'gandalf') library(nimbleHMC, lib.loc = '~/Documents/') else library(nimbleHMC)

