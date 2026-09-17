#!/bin/bash
#
#SBATCH -J loop.%j # le nom de votre job
#SBATCH -p long # la partition
#SBATCH -n 8 # nombre de processeurs
#SBATCH -N 1 # nombre de noeuds, DO NOT CHANGE THIS SETTING
#SBATCH --mem 60000 # mémoire (en MB)
#SBATCH --array=0
#SBATCH -o ./log/loop.%j.o # la sortie standard
#SBATCH -e ./log/loop.%j.e # la sortie erreur standard
#SBATCH --mail-type=END # notification par e-mail à la fin en cas de réussite ou erreur
#SBATCH --mail-user=xing.ma@universite-paris-saclay.fr # l’adresse e-mail

export NB_CORE=$SLURM_CPUS_ON_NODE # nb de processeurs alloués sur le noeud (spécifié par le parametre -n)
export OMP_NUM_THREADS=$NB_CORE  # limite du nombre thread à l'exécution
export OMP_THREAD_LIMIT=$NB_CORE # limite du nombre de thread pour l'ensemble du programme


mysamplePaths=($(echo *tag ))
f=${mysamplePaths[$SLURM_ARRAY_TASK_ID]}

cmd="makeTagDirectory $f -update -removeSelfLigation -removePEbg -genome $reference_genome"
echo $cmd
eval $cmd

cmd="analyzeHiC $f -cpu 8 -res 1000 -bgonly -norm"
echo $cmd
eval $cmd


cmd1="analyzeHiC $f -res 1000 -superRes 1000 -zscore 1.0 -nomatrix -cpu 8 -interactions `basename $f -tag`_1K_Z1.txt"
echo $cmd1
eval $cmd1

