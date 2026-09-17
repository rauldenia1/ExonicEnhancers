#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=20G
#SBATCH -t 01:00:00
#SBATCH -J tom_loopmaxZ
#SBATCH -o /NetScratch/CHROCYC/xma4/lab_project/Arab_microC/tomato/EE_ATAConly/loop/log/mz.%j.out
#SBATCH -e /NetScratch/CHROCYC/xma4/lab_project/Arab_microC/tomato/EE_ATAConly/loop/log/mz.%j.err
set -euo pipefail
module load R/4.1.0

DIR=/NetScratch/CHROCYC/xma4/lab_project/Arab_microC/tomato/EE_ATAConly/loop
cd "$DIR"

# highest Z-score of loops anchored by each gene
awk -F'\t' '{if($2>m[$1]) m[$1]=$2} END{for(g in m) printf "%s\t%.4f\n", g, m[g]}' gene_Z.tsv > gene_maxZ.txt
echo "genes with max loop Z: $(wc -l < gene_maxZ.txt)"

cat > tomato_loop_maxzscore_plot.R <<'RSCRIPT'
suppressMessages({library(dplyr); library(ggplot2)})
setwd("/NetScratch/CHROCYC/xma4/lab_project/Arab_microC/tomato/EE_ATAConly/loop")

fpkm <- read.table("fpkm_average_wt_order", header=TRUE, sep="\t", stringsAsFactors=FALSE)
colnames(fpkm) <- c("gene_id","fpkm")
qmap <- setNames(ntile(fpkm$fpkm, 5), fpkm$gene_id)
ee   <- scan("tom_EEgenes.txt", what="", quiet=TRUE)

z <- read.table("gene_maxZ.txt", sep="\t", col.names=c("gene_id","maxZ"), stringsAsFactors=FALSE)
z$grp <- ifelse(z$gene_id %in% ee, "EE gene",
         ifelse(z$gene_id %in% names(qmap), paste0("Q", qmap[z$gene_id]), NA))
z <- z[!is.na(z$grp), ]
disp <- c("Q1"="Q1\n(low)","Q2"="Q2","Q3"="Q3","Q4"="Q4","Q5"="Q5\n(high)","EE gene"="EE gene")
lev  <- c("Q1\n(low)","Q2","Q3","Q4","Q5\n(high)","EE gene")
z$group <- factor(disp[z$grp], levels=lev)

gcol <- c("Q1\n(low)"="#B0BEC5","Q2"="#90A4AE","Q3"="#78909C","Q4"="#607D8B","Q5\n(high)"="#455A64","EE gene"="#D4736A")
ymax <- max(z$maxZ); ytxt <- ymax*1.15
ee_vals <- z$maxZ[z$group=="EE gene"]; qg <- lev[1:5]
pv <- data.frame(group=factor(qg, levels=lev),
                 p=sapply(qg, function(g) wilcox.test(ee_vals, z$maxZ[z$group==g])$p.value))
pv$label <- ifelse(pv$p<2.2e-16,"p < 2.2e-16", paste0("p = ", formatC(pv$p, format="e", digits=2)))
pv$y <- ytxt

p <- ggplot(z, aes(group, maxZ, fill=group)) +
  geom_boxplot(outlier.shape=16, outlier.size=0.35, outlier.alpha=0.2, outlier.colour="grey50",
               width=0.6, linewidth=0.4) +
  scale_fill_manual(values=gcol) +
  scale_y_log10(name="Max Z-score of anchored loops (per gene)", limits=c(min(z$maxZ), ymax*1.35)) +
  geom_text(data=pv, aes(group, y, label=label), size=2.5, colour="grey30", inherit.aes=FALSE) +
  scale_x_discrete(name=NULL) + theme_bw() +
  theme(panel.grid.major.x=element_blank(), panel.grid.minor=element_blank(),
        panel.grid.major.y=element_line(colour="grey90", linewidth=0.3),
        legend.position="none", axis.text=element_text(size=10, colour="black"))
ggsave("fig_loop_maxzscore_by_group_tomato.pdf", p, width=5.5, height=4.5)
ggsave("fig_loop_maxzscore_by_group_tomato.png", p, width=5.5, height=4.5, dpi=150)

z %>% group_by(group) %>% summarise(n=n(), median_maxZ=round(median(maxZ),1)) %>% as.data.frame() %>% print()
RSCRIPT
Rscript tomato_loop_maxzscore_plot.R
echo DONE
