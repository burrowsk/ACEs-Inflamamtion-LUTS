#!/bin/bash

# 1 create a new dir and copy files to it
# perform the below
# perform the appending of the files to get 1 log
# reun twice once for each model (change postfix=".log" - postfix="_model2.log")

cd /user/work/epwkb/inflam

outcome=("UI" "Daywetting" "Bedwetting" "Frequency" "Nocturia" "Urgency" "Voiding_postponement" "Voiding_volume") 
exposure=("ACE_score")

for i in "${exposure[@]}"; do
	for j in "${outcome[@]}"; do

cd /user/home/epwkb/scratch/inflam/${i}/${j}/mod2

echo ${i}
echo ${j}

# create dir for new file names
mkdir filenames
cp *log filenames/

cd filenames/

# create dir for new filenames
mkdir filesuse

# This works https://stackoverflow.com/questions/55754/how-to-zero-pad-numbers-in-file-names-in-bash
prefix="g-form_${i}_${j}_"
postfix="_model2.log"
targetDir="filesuse"
paddingLength=2

for file in ${prefix}[0-9]*${postfix}; do
  # strip the prefix off the file name
  postfile=${file#$prefix}
  # strip the postfix off the file name
  number=${postfile%$postfix}
  # subtract 1 from the resulting number
  #i=$((number-1))
  # copy to a new name with padded zeros in a new folder
  cp ${file} "$targetDir"/$(printf $prefix%0${paddingLength}d$postfix $number)

done
done
done
