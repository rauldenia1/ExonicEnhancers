library(tidyverse)
library(GenomicRanges)
library(GenomicFeatures)
library(txdbmaker)
library(data.table)

# -------------------------------------------------------------------------
# 0. Global Setup & Plot Themes
# -------------------------------------------------------------------------
sample_fill_colors <- c("Permissive" = "#238b45", "Stringent" = "#2b8cbe")
sample_edge_colors <- c("Permissive" = "#00441b", "Stringent" = "#08589e")

occupancy_colors <- c(
  "0-25%"   = "#8B8CB3",
  "25-50%"  = "#8F737D",
  "50-75%"  = "#998665",
  "75-100%" = "#A2A79E"
)

theme_relative_density <- function() {
  theme_classic() +
    theme(
      plot.title      = element_text(face = "bold", size = 14),
      axis.text       = element_text(size = 10),
      legend.position = "top"
    )
}

scale_x_tss_relative <- function() {
  scale_x_continuous(
    breaks = c(0, 25, 50, 75, 100),
    labels = c("0%\n(TSS)", "25%", "50%", "75%", "100%")
  )
}

# -------------------------------------------------------------------------
# 1. Genome Annotation Setup
# -------------------------------------------------------------------------
gff_path <- "Tomato genome/SollycM82_genes_v1.1.1.gff3"

txdb <- makeTxDbFromGFF(
  file     = gff_path,
  format   = "auto",
  organism = "Solanum lycopersicum"
)

# Genome CDS (exons) annotation
genome_cds <- cds(txdb, columns = "gene_id")

# Extraer el ID primario si un CDS pertenece a múltiples genes
cds_gene_ids <- vapply(
  genome_cds$gene_id,
  function(x) if (length(x) > 0) x[1] else NA_character_,
  FUN.VALUE = character(1)
)
genome_cds$gene_id <- sub("^gene:", "", cds_gene_ids)
genome_cds$gene_id <- sub("_[0-9]+$", "", genome_cds$gene_id)

# Genome gene definitions
genome_genes <- genes(txdb)
genome_genes$gene_id <- sub("^gene:", "", as.character(genome_genes$gene_id))
genome_genes$gene_id <- sub("_[0-9]+$", "", genome_genes$gene_id)

genome_genes_df <- as.data.frame(genome_genes)

# -------------------------------------------------------------------------
# 2. Genomic Processing Function
# -------------------------------------------------------------------------
process_genomic_bed <- function(bed_path, sample_name) {
  # Read BED file and convert to GRanges
  bed_df <- fread(bed_path)
  colnames(bed_df)[1:3] <- c("seqnames", "start", "end")
  
  enhancers_gr <- makeGRangesFromDataFrame(
    bed_df,
    keep.extra.columns      = TRUE,
    starts.in.df.are.0based = TRUE
  )
  
  # Overlap enhancers with coding regions (CDS)
  hits <- findOverlaps(enhancers_gr, genome_cds, type = "within")
  
  matched_enhancers_df <- as.data.frame(enhancers_gr[queryHits(hits)])
  matched_cds_df       <- as.data.frame(genome_cds[subjectHits(hits)])
  matched_enhancers_df$true_gene_id <- matched_cds_df$gene_id
  
  # Extract midpoints and compute TSS distance
  enhancer_midpoints <- matched_enhancers_df$start + (matched_enhancers_df$width / 2)
  
  gene_indices <- match(matched_enhancers_df$true_gene_id, genome_genes_df$gene_id)
  gene_strands <- genome_genes_df$strand[gene_indices]
  gene_starts  <- genome_genes_df$start[gene_indices]
  gene_ends    <- genome_genes_df$end[gene_indices]
  gene_lengths <- genome_genes_df$width[gene_indices]
  
  matched_enhancers_df$distance_to_tss <- ifelse(
    gene_strands == "+",
    enhancer_midpoints - gene_starts,
    gene_ends - enhancer_midpoints
  )
  matched_enhancers_df$gene_length <- gene_lengths
  
  # Deduplicate by coordinates and calculate relative positions
  matched_enhancers_df |>
    distinct(seqnames, start, end, .keep_all = TRUE) |>
    mutate(
      relative_position     = distance_to_tss / gene_length,
      relative_position_pct = relative_position * 100,
      gene_coverage_pct     = (width / gene_length) * 100,
      occupancy = cut(
        gene_coverage_pct,
        breaks         = c(0, 25, 50, 75, 100),
        labels         = c("0-25%", "25-50%", "50-75%", "75-100%"),
        include.lowest = TRUE,
        right          = TRUE
      ),
      sample = sample_name
    )
}

# -------------------------------------------------------------------------
# 3. Load & Process Datasets
# -------------------------------------------------------------------------
files_to_process <- list(
  list(path = "tomato_EE_list/tomato_permissive_EEs.txt", sample = "Permissive"),
  list(path = "tomato_EE_list/tomato_stringent_EEs.txt",  sample = "Stringent")
)

all_ee_data <- files_to_process |>
  map_dfr(\(x) process_genomic_bed(x$path, x$sample))

# -------------------------------------------------------------------------
# 4. Plot 1: Relative Distance to TSS (All EEs)
# -------------------------------------------------------------------------
plot_ee_position <- ggplot(
  all_ee_data,
  aes(x = relative_position_pct, fill = sample, color = sample)
) +
  geom_density(alpha = 0.4, linewidth = 1) +
  scale_fill_manual(values = sample_fill_colors) +
  scale_color_manual(values = sample_edge_colors) +
  scale_x_tss_relative() +
  labs(
    x = "Relative position (%)",
    y = "EE density"
  ) +
  theme_relative_density()

print(plot_ee_position)
ggsave("ee_relative_position_comparison.svg", plot = plot_ee_position, width = 6, height = 5)

# -------------------------------------------------------------------------
# 5. Plot 2: Gene Occupancy / Coverage Stacked Barplot
# -------------------------------------------------------------------------
plot_ee_occupancy <- ggplot(all_ee_data, aes(x = sample, fill = occupancy)) +
  geom_bar(color = NA, alpha = 0.95, width = 0.5, position = "stack") +
  geom_text(
    stat  = "count",
    aes(label = after_stat(count)),
    vjust = -0.5,
    size  = 4
  ) +
  scale_fill_manual(values = occupancy_colors, drop = FALSE) +
  scale_y_continuous(
    limits = c(0, 1000),
    breaks = seq(0, 1000, by = 100),
    expand = expansion(mult = c(0, 0.1))
  ) +
  theme_minimal() +
  theme(
    panel.grid      = element_blank(),
    plot.title      = element_text(face = "bold", size = 13),
    axis.line       = element_line(color = "gray"),
    legend.position = "right"
  ) +
  labs(
    x    = NULL,
    y    = "Number of EEs",
    fill = "Occupancy"
  )

print(plot_ee_occupancy)
ggsave("ee_gene_occupancy_comparison.svg", plot = plot_ee_occupancy, width = 3.5, height = 5)

# -------------------------------------------------------------------------
# 6. Plot 3: Relative Distance (Coverage <= 25%)
# -------------------------------------------------------------------------
filtered_ee_data_25 <- all_ee_data |>
  filter(gene_coverage_pct <= 25)

plot_ee_position_25 <- ggplot(
  filtered_ee_data_25,
  aes(x = relative_position_pct, fill = sample, color = sample)
) +
  geom_density(alpha = 0.4, linewidth = 1) +
  scale_fill_manual(values = sample_fill_colors) +
  scale_color_manual(values = sample_edge_colors) +
  scale_x_tss_relative() +
  labs(
    x     = "Relative position (%)",
    y     = "EE density",
    title = "Coverage <= 25%"
  ) +
  theme_relative_density()

print(plot_ee_position_25)
ggsave("ee_relative_position_coverage_le25.svg", plot = plot_ee_position_25, width = 6, height = 5)

# -------------------------------------------------------------------------
# 7. Plot 4: Global Comparison (All vs <= 25% Coverage)
# -------------------------------------------------------------------------
combined_subset_data <- bind_rows(
  all_ee_data |> mutate(dataset_type = "All"),
  filtered_ee_data_25 |> mutate(dataset_type = "0-25%")
)

plot_all_vs_25 <- ggplot(
  combined_subset_data,
  aes(
    x        = relative_position_pct,
    color    = sample,
    linetype = dataset_type,
    group    = interaction(sample, dataset_type)
  )
) +
  geom_density(linewidth = 1, alpha = 0.15) +
  scale_color_manual(
    name   = "Condition",
    values = c("Permissive" = "#AACFE8", "Stringent" = "#636E8F")
  ) +
  scale_linetype_manual(
    name   = "Dataset",
    values = c("0-25%" = "solid", "All" = "dashed")
  ) +
  scale_x_tss_relative() +
  labs(
    x = "Relative position (%)",
    y = "EE density"
  ) +
  theme_classic() +
  theme(
    plot.title   = element_text(face = "bold", size = 14),
    axis.text    = element_text(size = 10),
    legend.position = "right",
    legend.box   = "vertical",
    legend.title = element_text(face = "bold", size = 9),
    legend.text  = element_text(size = 8)
  )

print(plot_all_vs_25)
ggsave("ee_all_vs_coverage_le25_comparison.svg", plot = plot_all_vs_25, width = 4.5, height = 4.5)