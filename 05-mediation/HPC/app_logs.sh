#!/bin/bash

# Append the mediation log files together for R
cd /user/work/epwkb/inflam

outcome=("UI" "Daywetting" "Bedwetting" "Frequency" "Nocturia" "Urgency" "Voiding_postponement" "Voiding_volume") 
exposure=("ACE_score")

for i in "${exposure[@]}"; do
	for j in "${outcome[@]}"; do

cd /user/home/epwkb/scratch/inflam/${i}/${j}/

echo ${i}
echo ${j}
grep "closed on:" *.log | wc -l

rm all_logs_${i}_${j}.txt 
cat *.log >> all_logs_${i}_${j}.txt
mv all_logs_${i}_${j}.txt ../../logs/

done
done