#!/bin/sh

clear

if [ -z "$1" ]
then
    # argument 1 is the empty string
    git diff
else
    # argument 1 is provided
    git diff $(git merge-base --fork-point "$1")
fi

