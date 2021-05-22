cd ~/github/nimble/nimble/packages/

rm -f nimble/inst/CppCode/*.o
rm -f nimble/inst/CppCode/*.so
rm -f nimble/src/*.o
rm -f nimble/src/*.so
rm -f nimble_*.tar.gz
 
R CMD build \
  --no-build-vignettes \
  --no-manual \
  --no-resave-data \
  nimble

if [ -z "$1" ]
then
    # argument 1 is the empty string
    R CMD INSTALL \
      --no-docs \
      --no-html \
      --no-data \
      --no-demo \
      --no-multiarch \
      --no-byte-compile \
      nimble_*.tar.gz
else
    # argument 1 is provided
    rm -rf ~/temp/builds/"$1"
    mkdir ~/temp/builds/"$1"
    R CMD INSTALL \
      --no-docs \
      --no-html \
      --no-data \
      --no-demo \
      --no-multiarch \
      --no-byte-compile \
      -l ~/temp/builds/"$1" \
      nimble_*.tar.gz
    ####VER=$(Rscript -e "cat(gsub('/Resources/.*$', '', gsub('^.*Versions/', '', .libPaths())))")
    ####cp  -r ~/temp/builds/"$1"/nimble /Library/Frameworks/R.framework/Versions/Current/Resources/library
    ####cp  -r ~/temp/builds/"$1"/nimble /Library/Frameworks/R.framework/Versions/"$VER"/Resources/library
fi




