library(tidyverse)
library(eulerr)
library(gplots)
library(writexl)
library(paletteer)
library(svglite)

# -------------------------------------------------------------------------
# 1. Orthogroup Identification & Matching
# -------------------------------------------------------------------------
# Load target gene lists
ee_genes_arabidopsis <- read_lines("geneIDs_EEs_Ara.txt")
ee_genes_tomato      <- read_lines("gene_IDs_permissive.txt")

# Build regex pattern for fast vectorized matching
ara_regex_pattern    <- paste0("\\b(", paste(ee_genes_arabidopsis, collapse = "|"), ")\\b")
tomato_regex_pattern <- paste0("\\b(", paste(ee_genes_tomato, collapse = "|"), ")\\b")

# Process OrthoFinder table
orthogroups_raw <- read_tsv("Orthogroups.tsv", show_col_types = FALSE)

orthogroups <- orthogroups_raw |>
  rename(
    orthogroup_id = 1,
    arabidopsis   = 2,
    tomato        = 3
  ) |>
  filter(
    !is.na(arabidopsis) & arabidopsis != "",
    !is.na(tomato) & tomato != ""
  ) |>
  mutate(
    # Normalize tomato gene transcript IDs (e.g., .1.1 to .1)
    tomato          = str_replace(tomato, "(\\.[0-9]+)\\.[0-9]+", "\\1"),
    has_ee_ara     = str_detect(arabidopsis, ara_regex_pattern),
    has_ee_tomato  = str_detect(tomato, tomato_regex_pattern),
    has_both       = has_ee_ara & has_ee_tomato
  )

write_csv2(orthogroups, "EEs_Orthogroups_permissive.csv")

# -------------------------------------------------------------------------
# 2. Venn Diagram & Intersection Export
# -------------------------------------------------------------------------
permissive_orthogroups <- read_csv2("EEs_Orthogroups_permissive.csv")

venn_input_list <- list(
  EEs_arabidopsis = permissive_orthogroups$Orthogroups[permissive_orthogroups$EEs_arabidopsis] |> na.omit(),
  EEs_tomato      = permissive_orthogroups$Orthogroups[permissive_orthogroups$ATAC_tomato] |> na.omit()
)

# Plot Venn Diagram
venn_fit <- euler(venn_input_list)

venn_plot <- plot(
  venn_fit,
  quantities = list(type = "counts", cex = 1.5, font = 2, col = "gray30"),
  shapes     = "circle",
  fills      = "white",
  edges      = list(col = c("#78C8F8FF", "#DD88ACFF", "#D0D0E0FF"), lex = 4),
  labels     = FALSE,
  legend     = list(side = "top", font = 1, cex = 1)
)

print(venn_plot)
ggsave("venn_arabidopsis_tomato_permissive.svg", plot = venn_plot, width = 5, height = 5)

# Export raw Venn intersection groups to Excel
venn_intersections <- attr(venn(venn_input_list, show.plot = FALSE), "intersections")
max_intersection_len <- max(lengths(venn_intersections))

intersections_df <- venn_intersections |>
  map(\(x) c(x, rep(NA, max_intersection_len - length(x)))) |>
  as_tibble()

write_xlsx(intersections_df, path = "venn_intersections_permissive.xlsx")

# -------------------------------------------------------------------------
# 3. Extract Common Genes from Shared Orthogroups
# -------------------------------------------------------------------------
permissive_ara_genes    <- read_lines("geneIDs_EEs_Ara.txt")
permissive_tomato_genes <- read_lines("gene_IDs_permissive.txt")

shared_orthogroups <- permissive_orthogroups |>
  filter(Both)

# Parse and filter Arabidopsis IDs
shared_ara_genes <- shared_orthogroups$Arabidopsis |>
  str_split(",\\s*") |>
  unlist() |>
  str_remove("\\..*") |>
  unique() |>
  intersect(permissive_ara_genes)

# Parse and filter Tomato IDs
shared_tomato_genes <- shared_orthogroups$Tomato |>
  str_split(",\\s*") |>
  unlist() |>
  unique() |>
  intersect(permissive_tomato_genes)

# Export paired table padded with NAs
max_gene_len <- max(length(shared_ara_genes), length(shared_tomato_genes))

common_genes_df <- tibble(
  Arabidopsis = `length<-`(shared_ara_genes, max_gene_len),
  Tomato      = `length<-`(shared_tomato_genes, max_gene_len)
)

write_csv(common_genes_df, "common_genes_permissive.tsv")

# -------------------------------------------------------------------------
# 4. Gene Ontology (GO) Enrichment Dotplot
# -------------------------------------------------------------------------
total_input_genes <- 190

go_plot_data <- read_tsv("GO_enrichment_common_permissive_190.txt", show_col_types = FALSE) |>
  mutate(GeneRatio = Count / total_input_genes) |>
  filter(Count < 100, `q-value` < 0.05) |>
  slice_min(order_by = `p-value`, n = 5, with_ties = FALSE)

min_ratio <- min(go_plot_data$GeneRatio)
max_ratio <- max(go_plot_data$GeneRatio)

go_enrichment_plot <- ggplot(go_plot_data, aes(x = GeneRatio, y = fct_reorder(Term, GeneRatio))) +
  geom_hline(yintercept = seq_len(nrow(go_plot_data)), color = "gray95", linewidth = 0.5) +
  geom_point(aes(size = Count, color = `p-value`)) +
  scale_color_paletteer_c("ggthemes::Blue Light", direction = -1) +
  scale_size_continuous(range = c(5, 15)) +
  scale_x_continuous(
    limits = c(min_ratio - 0.09, max_ratio + 0.09),
    breaks = seq(round(min_ratio, 2), round(max_ratio, 2), by = 0.1)
  ) +
  scale_y_discrete(labels = label_wrap_gen(25)) +
  theme_bw() +
  labs(
    x     = "GeneRatio",
    y     = NULL,
    color = "p-value",
    size  = "Count"
  ) +
  theme(
    panel.grid.major = element_line(color = "gray96"),
    panel.grid.minor = element_blank(),
    axis.text.y      = element_text(size = 12, color = "black"),
    axis.text.x      = element_text(size = 12, color = "black"),
    legend.text      = element_text(size = 8),
    legend.title     = element_text(size = 8, face = "bold"),
    aspect.ratio     = 1.5
  ) +
  guides(
    size = guide_legend(
      override.aes = list(
        shape  = 21,
        stroke = 0.5,
        color  = "black",
        fill   = "white"
      )
    )
  )

print(go_enrichment_plot)
ggsave(
  filename = "go_enrichment_common_permissive.svg",
  plot     = go_enrichment_plot,
  device   = svglite,
  width    = 10,
  height   = 5
)
