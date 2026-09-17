#!/bin/bash
#
#SBATCH -J homer.%j # le nom de votre job
#SBATCH -p chrocyc # la partition
#SBATCH -n 1 # nombre de processeurs
#SBATCH -N 1 # nombre de noeuds, DO NOT CHANGE THIS SETTING
#SBATCH --mem 30000 # mémoire (en MB)
#SBATCH --array=0
#SBATCH -o ./log/homer.%j.o # la sortie standard
#SBATCH -e ./log/homer.%j.e # la sortie erreur standard

export NB_CORE=$SLURM_CPUS_ON_NODE # nb de processeurs alloués sur le noeud (spécifié par le parametre -n)
export OMP_NUM_THREADS=$NB_CORE  # limite du nombre thread à l'exécution
export OMP_THREAD_LIMIT=$NB_CORE # limite du nombre de thread pour l'ensemble du programme

#awk '{print NR,$2,$3,$4,$5,$6,$7}' OFS="\t" /NetScratch/CHROCYC/xhe1/lab/2024/MicroC/HiC-pro/merge_rep/merged_results/HiC_result_tomato_3samples/data/Micro/Micro.allValidPairs > MicroC_allValidPairs.homer


mysamplePaths=($(echo *homer))
#mysamplePaths=($(echo ../merged_miC_Col.mapped.pairs))

file=${mysamplePaths[$SLURM_ARRAY_TASK_ID]}


#grep -v "#" $file| awk '{print NR,$2,$3,$6,$4,$5,$7}' OFS="\t" > `basename $file .mapped.pairs`.homer

#cmd0="awk 'BEGIN{FS=OFS=\"\t\"}{print NR,\$2,\$3,\$4,\$5,\$6,\$7}' $file > ./`basename $file .allValidPairs` _merged.homer"
#echo $cmd0
#eval $cmd0

cmd="makeTagDirectory $file-tag -format HiCsummary $file"
echo $cmd
eval $cmd

