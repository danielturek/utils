cd ~/github/nimble/nimbleHMC/
git pull
rm -f nimbleHMC_*.tar.gz
R CMD build nimbleHMC

rm -rf Documents/nimble

## now, in github/nimble/nimble repo,
## checkout branch ADoak_without_HMC, and run:
## build ADoak_without_HMC

cp -r ~/temp/builds/ADoak_without_HMC/nimble Documents/

R CMD INSTALL nimbleHMC_*.tar.gz -l ~/Documents/

## use in R scripts:
## if(Sys.info()['nodename'] == 'gandalf') library(nimbleHMC, lib.loc = '~/Documents/') else library(nimbleHMC)

