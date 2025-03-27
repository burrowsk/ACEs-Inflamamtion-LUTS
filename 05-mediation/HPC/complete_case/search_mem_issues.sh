#!/bin/bash

outcome=("UI" "Daywetting" "Bedwetting" "Frequency" "Nocturia" "Urgency" "Voiding_postponement" "Voiding_volume")
exposure=("ACE_score")

for i in "${exposure[@]}"; do for j in "${outcome[@]}"; do cd /user/home/epwkb/scratch/inflam/${i}/${j}/; echo ${i}; echo ${j}; rm mem_probs; grep ".unable to restore data" *.log >> /user/home/epwkb/scratch/inflam/mem_probs.txt;  done; done
