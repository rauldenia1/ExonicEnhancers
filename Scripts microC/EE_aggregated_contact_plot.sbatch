#!/bin/bash
#
#SBATCH -J hicAggregate # le nom de votre job
#SBATCH -p long # la partition
#SBATCH -n 2 # nombre de processeurs
#SBATCH -N 1 # nombre de noeuds, DO NOT CHANGE THIS SETTING
#SBATCH --mem 20000 # mémoire (en MB)
#SBATCH --array=0-4
#SBATCH -o ./log/agg.%j.o # la sortie standard
#SBATCH -e ./log/agg.%j.e # la sortie erreur standard
#SBATCH --mail-type=END # notification par e-mail à la fin en cas de réussite ou erreur
#SBATCH --mail-user=xing.ma@universite-paris-saclay.fr  # l’adresse e-mail

module load conda
conda activate hicexplorer

#Arabidopsis EE gene
hicAggregateContacts -m merged_miC_Col_5000.cool -out merged_miC_Col_5000_EE_gene.png --BED EE_gene.bed --mode intra-chr --transform  obs/exp --operationType mean --numberOfBins 16 --dpi 600

hicAggregateContacts -m merged_miC_Col_5000.cool -out merged_miC_Col_5000_EE_gene.png --BED EE_gene.bed --BED2 shuf_1k_exl_EE.bed --mode intra-chr --transform  obs/exp --operationType mean --numberOfBins 16 --dpi 600 --vMin 1.11 --vMax 1.28

#Tomato EE gene
hicAggregateContacts -m merged_Wt_5000.cool -out merged_Wt_5000_tomato_EE_gene.png --BED EE_gene_tomato.bed --mode intra-chr --transform  obs/exp --operationType mean --numberOfBins 16 --dpi 600 --vMin 1.5 --vMax 2.15

hicAggregateContacts -m merged_Wt_5000.cool -out merged_Wt_5000_tomato_EE_gene_vs_ctrl.png --BED EE_gene_tomato.bed --BED2 random_ctrl_tomato.bed --mode intra-chr --transform  obs/exp --operationType mean --numberOfBins 16 --dpi 600 --vMin 1.5 --vMax 2.15
