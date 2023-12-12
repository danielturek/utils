#!/bin/bash

## this open command needs to be first:
open -a "Google Chrome" http://localhost:3932/

Rscript -e "
library(blogdown)
siteDir <- '~/github/ds201/site/'
blogdown::serve_site(.site_dir = siteDir, port = 3932)
q('no')
"


