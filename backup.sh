#!/bin/bash

rm -rf ~/github/private/backup/*

cp -r ~/utils            ~/github/private/backup/
cp    ~/.bash_profile    ~/github/private/backup/
cp    ~/.emacs           ~/github/private/backup/
cp    ~/.gitconfig       ~/github/private/backup/
cp    ~/.rprofile        ~/github/private/backup/
cp    ~/scratch.R        ~/github/private/backup/

cd ~/github/private
~/utils/commit.sh automated file backup

