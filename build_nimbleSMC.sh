cd ~/github/nimble/nimbleSMC/packages/

rm -f nimbleSMC_*.tar.gz

R CMD build nimbleSMC

R CMD INSTALL nimbleSMC_*.tar.gz


