#!/bin/bash
#SBATCH -p chrocyc
#SBATCH -c 4
#SBATCH --mem=20G
#SBATCH -t 02:00:00
#SBATCH -J tom_strict
#SBATCH -o /NetScratch/CHROCYC/xma4/lab_project/Arab_microC/tomato/EE_ATAConly/loop/log/s.%j.out
#SBATCH -e /NetScratch/CHROCYC/xma4/lab_project/Arab_microC/tomato/EE_ATAConly/loop/log/s.%j.err
set -euo pipefail
module load bedtools/2.31.1
module load R/4.1.0

DIR=/NetScratch/CHROCYC/xma4/lab_project/Arab_microC/tomato/EE_ATAConly/loop
cd "$DIR"
LOOP=$DIR/WT_loops_intra_p0.01_z1.txt

# ---- strict loops: Z>1000 & FDR<0.05 & distance>10kb (Arabidopsis Z1k/>10kb definition) ----
awk -F'\t' 'NR==1 || ($17>1000 && $19<0.05 && $14>10000)' "$LOOP" > WT_loops_Z1k_10kb.txt
echo "strict loops (Z>1000,FDR<0.05,>10kb): $(( $(wc -l < WT_loops_Z1k_10kb.txt) - 1 ))"

# ---- anchors from strict loops ----
awk -F'\t' 'NR>1{print $3"\t"$4"\t"$5"\t"$1":1"; print $9"\t"$10"\t"$11"\t"$1":2"}' WT_loops_Z1k_10kb.txt \
  | sort -k1,1 -k2,2n > anchors_strict.bed
bedtools intersect -a anchors_strict.bed -b genes_sorted.bed -wa -wb | awk -F'\t' '{print $4"\t"$8}' > anchor_gene_strict.tsv

# ---- EE-anchored target genes (strict), excluding EE genes ----
bedtools intersect -a anchors_strict.bed -b ee_genes.bed -wa | awk -F'\t' '{print $4}' | sort -u > ee_anchor_tags_strict.txt
awk -F':' '{other=($2==1?2:1); print $1":"other}' ee_anchor_tags_strict.txt | sort -u > partner_tags_strict.txt
awk -F'\t' 'NR==FNR{p[$1]; next} ($1 in p){print $2}' partner_tags_strict.txt anchor_gene_strict.tsv | sort -u > ee_target_raw_strict.txt
comm -23 ee_target_raw_strict.txt ee_ids.txt > ee_target_genes_Z1k10kb.txt
echo "EE-anchored targets (strict, excl EE): $(wc -l < ee_target_genes_Z1k10kb.txt)"

# ---- expression boxplot: EE-anchored target vs non-anchored ----
cat > tomato_expr_strict.R <<'RSCRIPT'
suppressMessages({library(dplyr); library(ggplot2)})
setwd("/NetScratch/CHROCYC/xma4/lab_project/Arab_microC/tomato/EE_ATAConly/loop")
fpkm <- read.table("fpkm_average_wt_order", header=TRUE, sep="\t", stringsAsFactors=FALSE)
colnames(fpkm) <- c("gene_id","fpkm")
targets <- scan("ee_target_genes_Z1k10kb.txt", what="", quiet=TRUE)
fpkm$group <- ifelse(fpkm$gene_id %in% targets, "EE-anchored\ntarget (non-EE)", "Non-anchored")
fpkm$group <- factor(fpkm$group, levels=c("EE-anchored\ntarget (non-EE)","Non-anchored"))
wt <- wilcox.test(fpkm$fpkm[fpkm$group=="EE-anchored\ntarget (non-EE)"], fpkm$fpkm[fpkm$group=="Non-anchored"])
plab <- if (wt$p.value < 2.2e-16) "p < 2.2e-16" else paste0("p = ", formatC(wt$p.value, format="e", digits=2))
meds <- fpkm %>% group_by(group) %>% summarise(n=n(), median_fpkm=median(fpkm), .groups="drop")
p2 <- ggplot(fpkm, aes(group, fpkm+1, fill=group)) +
  geom_boxplot(outlier.shape=16, outlier.size=0.3, outlier.alpha=0.15, outlier.colour="grey50", width=0.55, linewidth=0.4) +
  scale_fill_manual(values=c("EE-anchored\ntarget (non-EE)"="#D4736A","Non-anchored"="#90A4AE")) +
  scale_y_log10(name="log10 FPKM+1") +
  annotate("text", x=1.5, y=max(fpkm$fpkm), label=plab, size=4, colour="grey30") +
  scale_x_discrete(name=NULL) + theme_bw() +
  theme(panel.grid.major.x=element_blank(), panel.grid.minor=element_blank(),
        panel.grid.major.y=element_line(colour="grey90", linewidth=0.3),
        legend.position="none", axis.text=element_text(size=10, colour="black"))
ggsave("fig_expr_EE_targets_tomato_Z1k10kb.pdf", p2, width=4, height=4.5)
ggsave("fig_expr_EE_targets_tomato_Z1k10kb.png", p2, width=4, height=4.5, dpi=150)
cat("=== strict EE-target expression ===\n"); print(as.data.frame(meds)); cat(sprintf("Wilcoxon p = %.3g\n", wt$p.value))
RSCRIPT
Rscript tomato_expr_strict.R
echo DONE
