#!/bin/bash

##rm -rf ~/github/private/backup/*

files=".bash_profile .emacs .gitconfig .rprofile scratch.R"

for file in $files; do
    cp ~/$file ~/github/utils/
done

cd ~/github/utils
~/github/utils/commit.sh automated file backup

for file in $files; do
    cp ~/github/utils/$file ~/
done




#cp -r ~/utils            ~/github/private/backup/
#cp    ~/.bash_profile    ~/github/private/backup/
#cp    ~/.emacs           ~/github/private/backup/
#cp    ~/.gitconfig       ~/github/private/backup/
#cp    ~/.rprofile        ~/github/private/backup/
#cp    ~/scratch.R        ~/github/private/backup/

