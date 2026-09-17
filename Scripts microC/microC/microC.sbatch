#!/bin/bash
#
#SBATCH -J micro-C_hic # le nom de votre job
#SBATCH -p chrocyc # la partition
#SBATCH -n 8 # nombre de processeurs
#SBATCH -N 1 # nombre de noeuds, DO NOT CHANGE THIS SETTING
#SBATCH --mem 50G
#SBATCH --array=0-1
#SBATCH -o ./log/micro-C_hic.%j.o # la sortie standard
#SBATCH -e ./log/micro-C_hic.%j.e # la sortie erreur standard
#SBATCH --mail-type=END # notification par e-mail à la fin en cas de réussite ou erreur
#SBATCH --mail-user=xing.ma@universite-paris-saclay.fr # l’adresse e-mail


module load conda
conda activate pairtools_env

module load samtools
module load bwa
module load java


mysamplePaths=($(echo rawdata/*_R1.fq.gz ))
file=${mysamplePaths[$SLURM_ARRAY_TASK_ID]}
outputDir=.


#Map data with bwa mem

cmd1="bwa mem -5SP -T0 -t8 /NetScratch/CHROCYC/xma4/genome/TAIR10/chr_TAIR10/TAIR10.fasta $file `dirname $file`/`basename $file _R1.fq.gz`_R2.fq.gz -o $outputDir/`basename $file _R1.fq.gz`.sam"
echo $cmd1
eval $cmd1

f2=$outputDir/`basename $file _R1.fq.gz`.sam


#pairtools parse, parse .sam/.bam into pairs

cmd2="pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in 8 --nproc-out 8 --chroms-path /NetScratch/CHROCYC/xma4/genome/TAIR10/chr_TAIR10/TAIR10.sizes.chr  $f2 > $outputDir/`basename $f2`_parsed.pairsam"
echo $cmd2
eval $cmd2


#pairtools sort
#Sorting makes it more organised, and also helps to perform deduplication which cannot work on non-sorted files.

mkdir -p $outputDir/temp

cmd3="pairtools sort --nproc 8 --tmpdir=$outputDir/temp/  `basename $f2`_parsed.pairsam > $outputDir/`basename $f2`_sorted.pairsam"
echo $cmd3
eval $cmd3

#pairtools dedup

cmd4="pairtools dedup --nproc-in 8 --nproc-out 8 --mark-dups --output-stats $outputDir/`basename $f2`_stats.txt --output $outputDir/`basename $f2`_dedup.pairsam $outputDir/`basename $f2`_sorted.pairsam"
echo $cmd4
eval $cmd4


#split pairsam into pairs and sam files, the pairs file is used for downstream analysis

cmd5="pairtools split --nproc-in 8 --nproc-out 8 --output-pairs $outputDir/`basename $f2 sam`mapped.pairs --output-sam $outputDir/`basename $f2 sam`mapped.pairs.sam $outputDir/`basename $f2`_dedup.pairsam "
echo $cmd5
eval $cmd5


#convert pairs file into .hic file, the .hic file is a binary file that contains the contact matrix and other information, it can be visualized by juicebox or other tools
#java is needed to run juicer_tools.jar

cmd6="java -Xmx50g -Djava.awt.headless=true -jar /NetScratch/CHROCYC/xhe1/tools/juicer/CPU/common/juicer_tools.jar pre --threads 8 `basename $f2 sam`mapped.pairs  `basename $f2 sam`hic /NetScratch/CHROCYC/xma4/genome/TAIR10/chr_TAIR10/TAIR10.sizes.chr"
echo $cmd6
eval $cmd6


