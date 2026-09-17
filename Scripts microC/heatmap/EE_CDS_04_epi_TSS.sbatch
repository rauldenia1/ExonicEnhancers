#!/bin/bash
#SBATCH -J EEcds_epiTSS
#SBATCH -p short
#SBATCH -n 8
#SBATCH -N 1
#SBATCH --mem 14000
#SBATCH -o /NetScratch/CHROCYC/xma4/lab_project/Arab_microC/tomato/EE_CDS_ATAConly/epi_profile_TSS/log/dpts.%j.o
#SBATCH -e /NetScratch/CHROCYC/xma4/lab_project/Arab_microC/tomato/EE_CDS_ATAConly/epi_profile_TSS/log/dpts.%j.e
# ==============================================================================

set -e
export NB_CORE=$SLURM_CPUS_ON_NODE
export OMP_NUM_THREADS=$NB_CORE
export OMP_THREAD_LIMIT=$NB_CORE

module load deeptools

computeMatrix reference-point \
     -S $BW1/ATAC_HS_0h.bigwig \
        $BW1/coverage_H3K4me1-Rep2.bigwig \
        $BW1/coverage_H3K4me2-Rep2.bigwig \
        $BW1/coverage_H3K4me3-Rep2.bigwig \
        $BW1/coverage_WT_K9ac.bigwig \
        $BW1/coverage_Rim_H3K18ac_R2_S10.bigwig \
        $BW1/coverage_H3K36me3-Rep2.bigwig \
        $BW1/coverage_WT_K27me3.bigwig \
        $BW1/coverage_WT_K9me2.bigwig \
     -R $REGION --referencePoint TSS \
     -a 2000 -b 2000 --missingDataAsZero \
     -o EE_TSS_epi_matrix.gz -p 8

LAB=(ATAC H3K4me1 H3K4me2 H3K4me3 H3K9ac H3K18ac H3K36me3 H3K27me3 H3K9me2)
N=${#LAB[@]}
HCOL2=(); PCOL2=(); ZMIN=()
for ((i=0;i<N;i++)); do HCOL2+=("white,#2166AC"); PCOL2+=("#2166AC"); ZMIN+=("0"); done

for ext in pdf png; do
  dpi=600; [ "$ext" = png ] && dpi=300
  plotHeatmap -m EE_TSS_epi_matrix.gz -o EE_TSS_epi_heatmap.$ext --dpi $dpi \
       --samplesLabel "${LAB[@]}" \
       --plotTitle '' --whatToShow "plot, heatmap and colorbar" \
       --colorList "${HCOL2[@]}" \
       --regionsLabel "" --refPointLabel 'TSS' --xAxisLabel '' \
       --zMin "${ZMIN[@]}" \
       --sortUsing mean --sortUsingSamples 1 \
       --boxAroundHeatmaps no \
       --heatmapHeight 15 --heatmapWidth 4

done
echo "set 2 done"
echo DONE
