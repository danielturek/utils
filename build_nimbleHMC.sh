cd ~/github/nimble/nimbleHMC/

rm -f nimbleHMC_*.tar.gz

R CMD build nimbleHMC

R CMD INSTALL nimbleHMC_*.tar.gz


