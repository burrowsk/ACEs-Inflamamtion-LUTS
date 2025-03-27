#!/bin/bash

cd /user/work/epwkb/inflam

outcome=("UI" "Daywetting" "Bedwetting" "Frequency" "Nocturia" "Urgency" "Voiding_postponement" "Voiding_volume") 
exposure=("ACE_score" "Physical_abuse" "Sexual_abuse" "Emotional_abuse" "Emotional_neglect" "Bullying" "Violence_parents" "Substance_abuse" "Mental_hlth_suicide" "Convicted_offences" "Parental_separation")

for i in "${exposure[@]}"; do
	for j in "${outcome[@]}"; do

cd /user/home/epwkb/scratch/inflam/${i}/${j}

echo ${i}
echo ${j}


#mkdir mod2
#mv *model2.log mod2/

#mkdir out
#mv *.out out/

grep "closed on:" filenames/filesuse/*.log | wc -l

#rm all_logs_${i}_${j}.txt 
cat filenames/filesuse/*.log >> filenames/filesuse/all_logs_${i}_${j}.txt
mv filenames/filesuse/all_logs_${i}_${j}.txt ../../logs

grep "closed on:" mod2/filenames/filesuse/*.log | wc -l

#rm all_logs_${i}_${j}_mod2.txt 
cat mod2/filenames/filesuse/*.log >> mod2/filenames/filesuse/all_logs_${i}_${j}_mod2.txt
mv mod2/filenames/filesuse/all_logs_${i}_${j}_mod2.txt ../../logs


done
done