#!/bin/bash
# trim
#
#SBATCH -J merge.%j # le nom de votre job
#SBATCH -p chrocyc # la partition
#SBATCH -n 8 # nombre de processeurs
#SBATCH -N 1 # nombre de noeuds, DO NOT CHANGE THIS SETTING
#SBATCH --mem 50000 # mémoire (en MB)
#SBATCH --array=0
#SBATCH -o ./log/merge.%j.o # la sortie standard
#SBATCH -e ./log/merge.%j.e # la sortie erreur standard
#SBATCH --mail-type=END # notification par e-mail à la fin en cas de réussite ou erreur
#SBATCH --mail-user=xing.ma@universite-paris-saclay.fr # l’adresse e-mail

module load conda
conda activate pairtools_env

#merge pairs files from different replicates into one file

pairtools merge --nproc 8 --memory 50G miC_Col_rep1.mapped.pairs miC_Col_rep2.mapped.pairs -o merged_miC_Col.mapped.pairs


#convert pairs file into .hic file

module load java

java -Xmx50g -Djava.awt.headless=true -jar /NetScratch/CHROCYC/xhe1/tools/juicer/CPU/common/juicer_tools.jar pre --threads 8 merged_miC_Col.mapped.pairs merged_miC_Col.res.hic /NetScratch/CHROCYC/xma4/genome/TAIR10/chr_TAIR10/TAIR10.sizes.chr -r 500,1000,2000,5000,10000,25000,50000,100000

