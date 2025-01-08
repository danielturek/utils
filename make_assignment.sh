#!/bin/bash

# Check if exactly one argument is provided
if [ $# -ne 1 ]; then
  echo "Usage: $0 X-Y"
  exit 1
fi

# Get the input argument
input=$1

# Split the argument into X and Y
IFS='-' read -r X Y <<< "$input"

# Set the folder variable based on the value of X
case "$X" in
  hw)
    folder="hw"
    ;;
  lab)
    folder="labs"
    ;;
  ex)
    folder="exercises"
    ;;
  *)
    echo "Invalid value for X. Expected 'hw', 'lab', or 'ex'."
    exit 1
    ;;
esac

# Print the values of folder and the original argument
##echo "folder: $folder"
##echo "original argument: $input"


dir=$(pwd)
cp -r ~/github/ds201/ds201/starters/"$folder"/"$input" ~/github/ds201/starters/"$folder"/
cd ~/github/ds201/starters/"$folder"/"$input"/
git init
git remote add origin https://github.com/DS201-S25/"$input".git
git add --all
git commit -m"added starter repo materials"
git push -u origin main
cd "$dir"





