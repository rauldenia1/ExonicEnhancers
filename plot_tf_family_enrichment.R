library(tidyverse)

# -------------------------------------------------------------------------
# 0. Global Setup & Plot Themes
# -------------------------------------------------------------------------
condition_colors <- c(
  "Permissive" = "#AACFE8",
  "Stringent"  = "#636E8F"
)

theme_tf_barplot <- function() {
  theme_minimal() +
    theme(
      panel.grid      = element_blank(),
      axis.line.x     = element_line(color = "gray"),
      axis.line.y     = element_line(color = "gray"),
      axis.text.y     = element_text(size = 8),
      axis.text.x     = element_text(size = 9, color = "black", face = "bold", angle = 90, hjust = 1, vjust = 0.5),
      legend.position = "top",
      legend.title    = element_blank()
    )
}

# -------------------------------------------------------------------------
# 1. Load Reference Metadata
# -------------------------------------------------------------------------
metadata <- read_tsv("Sly_TF_binding_motifs_information.txt", show_col_types = FALSE)

# Reference family total counts
family_reference_totals <- metadata |>
  drop_na(Family) |>
  count(Family, name = "total_in_system")

# -------------------------------------------------------------------------
# 2. Processing Pipeline
# -------------------------------------------------------------------------
process_tfbs <- function(file_path, condition_name) {
  read_tsv(file_path, show_col_types = FALSE, comment = "#") |>
    filter(QVALUE < 0.05) |>
    left_join(
      metadata |> select(Gene_id, Family),
      by = c("ID" = "Gene_id")
    ) |>
    drop_na(Family) |>
    count(Family, name = "tf_count") |>
    mutate(condition = condition_name)
}

enrichment_files <- list(
  list(path = "TFBS enrichment/TFBS_enrichment_permissive.tsv", condition = "Permissive"),
  list(path = "TFBS enrichment/TFBS_enrichment_stringent.tsv",  condition = "Stringent")
)

# Complete zeros, filter by size, and compute ratios
comparison_data <- enrichment_files |>
  map_dfr(\(x) process_tfbs(x$path, x$condition)) |>
  complete(Family, condition, fill = list(tf_count = 0)) |>
  left_join(family_reference_totals, by = "Family") |>
  filter(total_in_system > 1) |>
  mutate(
    success_percentage = (tf_count / total_in_system) * 100,
    ratio_label        = paste0(tf_count, "/", total_in_system),
    # Reorder Family factor globally by mean success percentage across conditions
    Family             = fct_reorder(Family, success_percentage, .fun = mean, .desc = TRUE)
  )

# -------------------------------------------------------------------------
# 3. Barplot Visualization
# -------------------------------------------------------------------------
tf_family_plot <- ggplot(
  comparison_data,
  aes(x = Family, y = success_percentage, fill = condition)
) +
  geom_col(
    position = position_dodge(width = 0.8),
    color    = "white",
    width    = 0.7
  ) +
  geom_text(
    aes(label = ratio_label),
    position = position_dodge(width = 0.8),
    vjust    = 0.25,
    size     = 2.5,
    color    = "black",
    angle    = 90,
    hjust    = -0.1
  ) +
  scale_fill_manual(values = condition_colors) +
  scale_y_continuous(
    limits = c(0, 130),
    breaks = seq(0, 100, by = 25),
    expand = expansion(mult = c(0, 0.05))
  ) +
  theme_tf_barplot() +
  labs(
    x = NULL,
    y = "% (TFBS / Total in family)"
  )

print(tf_family_plot)
ggsave("tf_family_enrichment_comparison.svg", plot = tf_family_plot, width = 4, height = 4.5)