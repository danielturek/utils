#!/bin/bash

## this open command needs to be first:
open -a "Google Chrome" http://localhost:1000/

Rscript -e "
(source('~/temp/shiny/app-3/app.R'))
q('no')
"


