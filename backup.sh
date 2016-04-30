#!/bin/bash

files=".bash_profile .emacs .gitconfig .rprofile scratch.R"


for file in $files; do
    cp ~/$file ~/github/utils/
done


cd ~/github/utils
~/github/utils/commit.sh automated file backup


for file in $files; do
    cp ~/github/utils/$file ~/
done

