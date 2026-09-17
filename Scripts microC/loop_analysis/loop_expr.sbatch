#!/bin/bash
#SBATCH -p chrocyc
#SBATCH -c 8
#SBATCH --mem=30G
#SBATCH -t 04:00:00
#SBATCH -J tom_loopexpr
#SBATCH -o /NetScratch/CHROCYC/xma4/lab_project/Arab_microC/tomato/EE_ATAConly/loop/log/j.%j.out
#SBATCH -e /NetScratch/CHROCYC/xma4/lab_project/Arab_microC/tomato/EE_ATAConly/loop/log/j.%j.err
set -euo pipefail
module load bedtools/2.31.1
module load R/4.1.0

DIR=/NetScratch/CHROCYC/xma4/lab_project/Arab_microC/tomato/EE_ATAConly/loop
cd "$DIR"; mkdir -p log
GENEBED=$DIR/M82_gene.bed7
FPKM=$DIR/fpkm_average_wt_order
EEGENES=$DIR/tom_EEgenes.txt
LOOP=$DIR/WT_loops_intra_p0.01_z1.txt

# ---- sanitize gene bed (chr1-12, numeric coords, strip CR) ----
awk -F'\t' '{sub(/\r$/,"")} $1~/^chr[0-9]+$/ && $2~/^[0-9]+$/ && $3~/^[0-9]+$/' "$GENEBED" | sort -k1,1 -k2,2n > genes_sorted.bed

# ---- anchors: two per loop, tagged loopID:end ----
awk -F'\t' 'NR>1{print $3"\t"$4"\t"$5"\t"$1":1"; print $9"\t"$10"\t"$11"\t"$1":2"}' "$LOOP" | sort -k1,1 -k2,2n > anchors.bed
echo "anchors: $(wc -l < anchors.bed)"

# ---- anchor -> gene (tag, gene) ----
bedtools intersect -a anchors.bed -b genes_sorted.bed -wa -wb | awk -F'\t' '{print $4"\t"$8}' > anchor_gene.tsv

# ---- loop count per gene (distinct loops anchored per gene) ----
awk -F'\t' '{split($1,b,":"); print $2"\t"b[1]}' anchor_gene.tsv | sort -u | cut -f1 | sort | uniq -c \
  | awk '{print $2"\t"$1}' > loop_count_per_gene.txt
echo "genes with >=1 loop: $(wc -l < loop_count_per_gene.txt)"

# ---- EE-anchored target genes (partner anchor of EE-gene anchors), excluding EE genes ----
awk -F'\t' 'NR==FNR{k[$1]; next} ($4 in k)' "$EEGENES" genes_sorted.bed > ee_genes.bed
bedtools intersect -a anchors.bed -b ee_genes.bed -wa | awk -F'\t' '{print $4}' | sort -u > ee_anchor_tags.txt
awk -F':' '{other=($2==1?2:1); print $1":"other}' ee_anchor_tags.txt | sort -u > partner_tags.txt
awk -F'\t' 'NR==FNR{p[$1]; next} ($1 in p){print $2}' partner_tags.txt anchor_gene.tsv | sort -u > ee_target_raw.txt
sort -u "$EEGENES" > ee_ids.txt
comm -23 ee_target_raw.txt ee_ids.txt > ee_target_genes.txt
echo "EE-anchored target genes (excl EE): $(wc -l < ee_target_genes.txt)"

# ================== PLOTS (R) ==================
cat > tomato_loop_expr_plots.R <<'RSCRIPT'
suppressMessages({library(dplyr); library(ggplot2)})
setwd("/NetScratch/CHROCYC/xma4/lab_project/Arab_microC/tomato/EE_ATAConly/loop")

fpkm <- read.table("fpkm_average_wt_order", header=TRUE, sep="\t", stringsAsFactors=FALSE)
colnames(fpkm) <- c("gene_id","fpkm")
fpkm$Quintile <- ntile(fpkm$fpkm, 5)
ee <- scan("tom_EEgenes.txt", what="", quiet=TRUE)
qmap <- setNames(fpkm$Quintile, fpkm$gene_id)

## ---------- Analysis 1: loops per gene across expression groups + EE ----------
counts <- read.table("loop_count_per_gene.txt", sep="\t",
                     col.names=c("gene_id","n_loops"), stringsAsFactors=FALSE)
counts$grp <- ifelse(counts$gene_id %in% ee, "EE gene",
              ifelse(counts$gene_id %in% names(qmap), paste0("Q", qmap[counts$gene_id]), NA))
counts <- counts[!is.na(counts$grp), ]
disp <- c("Q1"="Q1\n(low)","Q2"="Q2","Q3"="Q3","Q4"="Q4","Q5"="Q5\n(high)","EE gene"="EE gene")
lev  <- c("Q1\n(low)","Q2","Q3","Q4","Q5\n(high)","EE gene")
counts$group <- factor(disp[counts$grp], levels=lev)

gcol <- c("Q1\n(low)"="#B0BEC5","Q2"="#90A4AE","Q3"="#78909C","Q4"="#607D8B","Q5\n(high)"="#455A64","EE gene"="#D4736A")
ymax <- max(counts$n_loops); ytxt <- ymax*1.15
ee_vals <- counts$n_loops[counts$group=="EE gene"]
qg <- c("Q1\n(low)","Q2","Q3","Q4","Q5\n(high)")
pdf1 <- data.frame(group=factor(qg, levels=lev),
                   p=sapply(qg, function(g) wilcox.test(ee_vals, counts$n_loops[counts$group==g])$p.value))
pdf1$label <- ifelse(pdf1$p<2.2e-16,"p < 2.2e-16", paste0("p = ", formatC(pdf1$p, format="e", digits=2)))
pdf1$y <- ytxt

p1 <- ggplot(counts, aes(group, n_loops, fill=group)) +
  geom_boxplot(outlier.shape=16, outlier.size=0.35, outlier.alpha=0.2, outlier.colour="grey50",
               width=0.6, linewidth=0.4) +
  scale_fill_manual(values=gcol) +
  scale_y_log10(name="Loops anchored per gene", limits=c(1, ymax*1.35)) +
  geom_text(data=pdf1, aes(group, y, label=label), size=2.5, colour="grey30", inherit.aes=FALSE) +
  scale_x_discrete(name=NULL) + theme_bw() +
  theme(panel.grid.major.x=element_blank(), panel.grid.minor=element_blank(),
        panel.grid.major.y=element_line(colour="grey90", linewidth=0.3),
        legend.position="none", axis.text=element_text(size=10, colour="black"))
ggsave("fig_loop_count_by_group_tomato.pdf", p1, width=5.5, height=4.5)
ggsave("fig_loop_count_by_group_tomato.png", p1, width=5.5, height=4.5, dpi=150)

## ---------- Analysis 2: expression of EE-anchored targets vs non-anchored ----------
targets <- scan("ee_target_genes.txt", what="", quiet=TRUE)
fpkm$group <- ifelse(fpkm$gene_id %in% targets, "EE-anchored\ntarget (non-EE)", "Non-anchored")
fpkm$group <- factor(fpkm$group, levels=c("EE-anchored\ntarget (non-EE)","Non-anchored"))
wt <- wilcox.test(fpkm$fpkm[fpkm$group=="EE-anchored\ntarget (non-EE)"],
                  fpkm$fpkm[fpkm$group=="Non-anchored"])
plab <- if (wt$p.value < 2.2e-16) "p < 2.2e-16" else paste0("p = ", formatC(wt$p.value, format="e", digits=2))
meds <- fpkm %>% group_by(group) %>% summarise(n=n(), median_fpkm=median(fpkm), .groups="drop")

p2 <- ggplot(fpkm, aes(group, fpkm+1, fill=group)) +
  geom_boxplot(outlier.shape=16, outlier.size=0.3, outlier.alpha=0.15, outlier.colour="grey50",
               width=0.55, linewidth=0.4) +
  scale_fill_manual(values=c("EE-anchored\ntarget (non-EE)"="#D4736A","Non-anchored"="#90A4AE")) +
  scale_y_log10(name="log10 FPKM+1") +
  annotate("text", x=1.5, y=max(fpkm$fpkm), label=plab, size=4, colour="grey30") +
  scale_x_discrete(name=NULL) + theme_bw() +
  theme(panel.grid.major.x=element_blank(), panel.grid.minor=element_blank(),
        panel.grid.major.y=element_line(colour="grey90", linewidth=0.3),
        legend.position="none", axis.text=element_text(size=10, colour="black"))
ggsave("fig_expr_EE_targets_tomato.pdf", p2, width=4, height=4.5)
ggsave("fig_expr_EE_targets_tomato.png", p2, width=4, height=4.5, dpi=150)

cat("=== Analysis 1: loops per gene (n per group) ===\n"); print(table(counts$group))
cat("=== Analysis 2: expression groups ===\n"); print(as.data.frame(meds))
cat(sprintf("Analysis2 Wilcoxon p = %.3g\n", wt$p.value))
RSCRIPT
Rscript tomato_loop_expr_plots.R
echo "DONE"
