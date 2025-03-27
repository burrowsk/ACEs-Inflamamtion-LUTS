#!/bin/bash

array=$1
outcome=("UI" "Daywetting" "Bedwetting" "Frequency" "Nocturia" "Urgency" "Voiding_postponement" "Voiding_volume") 
exposure=("ACE_score")

for i in "${exposure[@]}"; do
	for j in "${outcome[@]}"; do
#mkdir /user/home/epwkb/scratch/inflam/complete_case/${i}/
#mkdir /user/home/epwkb/scratch/inflam/complete_case/${i}/${j}/

cd /user/home/epwkb/scratch/inflam/complete_case/${i}/${j}/
#rm *log *out *txt
file=/user/home/epwkb/scratch/inflam/complete_case/submit_imps.sh
sed -i 's/#SBATCH --array=1/#SBATCH --array='$array'/' $file
sbatch /user/home/epwkb/scratch/inflam/complete_case/submit_imps.sh ${i} ${j}

done
done

