cd ~/github/nimble/nimbleHMC/
git pull
rm -f nimbleHMC_*.tar.gz
R CMD build nimbleHMC

rm -rf Documents/nimble

cd ~/github/nimble/nimble/
git checkout ADoak_without_HMC
build ADoak_without_HMC

cp -r ~/temp/builds/ADoak_without_HMC/nimble Documents/

cd ~/github/nimble/nimbleHMC/
R CMD INSTALL nimbleHMC_*.tar.gz -l ~/Documents/

## use in R scripts:
## if(Sys.info()['nodename'] == 'gandalf') library(nimbleHMC, lib.loc = '~/Documents/') else library(nimbleHMC)

