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

R CMD INSTALL \
  --no-docs \
  --no-html \
  --no-data \
  --no-help \
  --no-demo \
  --no-multiarch \
  --no-byte-compile \
  nimble_*.tar.gz
