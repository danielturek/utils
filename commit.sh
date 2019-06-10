#!/bin/bash

git pull
git add --all :/

#git commit -m"$*."
if [ -z "$1" ]
then
    git commit -m"updates"
else
    git commit -m "$*"
fi

git push
