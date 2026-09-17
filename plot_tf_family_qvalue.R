library(tidyverse)

# -------------------------------------------------------------------------
# 0. Global Setup & Plot Theme
# -------------------------------------------------------------------------
bar_fill_color <- "#636E8F"

theme_tf_barplot <- function() {
  theme_minimal() +
    theme(
      panel.grid  = element_blank(),
      axis.line.x = element_line(color = "gray"),
      axis.line.y = element_line(color = "gray"),
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(size = 10, color = "black", face = "bold", angle = 90, hjust = 1, vjust = 0.5)
    )
}

# -------------------------------------------------------------------------
# 1. Load Data
# -------------------------------------------------------------------------
metadata <- read_tsv(
  "Sly_TF_binding_motifs_information.txt",
  show_col_types = FALSE
)

tfbs_enrichment <- read_tsv(
  "TFBS enrichment/TFBS_enrichment_stringent.tsv",
  show_col_types = FALSE
)

# -------------------------------------------------------------------------
# 2. Data Processing & Family Normalization
# -------------------------------------------------------------------------
# Total TF reference counts per family in the reference genome
family_reference_totals <- metadata |>
  drop_na(Family) |>
  count(Family, name = "total_in_system")

# Identify enriched TFs and join family information
plot_data <- tfbs_enrichment |>
  filter(QVALUE < 0.05) |>
  left_join(
    metadata |> select(Gene_id, Family),
    by = c("ID" = "Gene_id")
  ) |>
  drop_na(Family) |>
  count(Family, name = "detected_tfs") |>
  left_join(family_reference_totals, by = "Family") |>
  filter(total_in_system > 1) |>
  mutate(
    success_percentage = (detected_tfs / total_in_system) * 100,
    ratio_label        = paste0(detected_tfs, " / ", total_in_system),
    Family             = fct_reorder(Family, success_percentage, .desc = TRUE)
  )

# -------------------------------------------------------------------------
# 3. Barplot Visualization
# -------------------------------------------------------------------------
tf_enrichment_plot <- ggplot(plot_data, aes(x = Family, y = success_percentage)) +
  geom_col(
    fill  = bar_fill_color,
    color = "white",
    width = 0.9
  ) +
  geom_text(
    aes(label = ratio_label),
    size    = 3,
    color   = "black",
    nudge_y = 10,
    angle   = 90
  ) +
  scale_y_continuous(
    limits = c(0, 110),
    breaks = seq(0, 100, by = 25),
    expand = expansion(mult = c(0, 0.1))
  ) +
  theme_tf_barplot() +
  labs(
    x = NULL,
    y = "% (TFBS / Total in family)"
  )

print(tf_enrichment_plot)
ggsave("tfbs_stringent_qvalue_enrichment.svg", plot = tf_enrichment_plot, width = 2.5, height = 4)