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
