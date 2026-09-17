library(tidyverse)

# -------------------------------------------------------------------------
# 0. Global Setup & Plot Themes
# -------------------------------------------------------------------------
condition_colors <- c(
  "Permissive" = "#AACFE8",
  "Stringent"  = "#636E8F"
)

# Custom minimal theme to avoid code duplication across plots
theme_ee_distribution <- function() {
  theme_minimal() +
    theme(
      panel.grid   = element_blank(),
      axis.line.x  = element_line(color = "gray"),
      axis.line.y  = element_line(color = "gray"),
      axis.text.y  = element_text(size = 8),
      axis.text.x  = element_text(size = 8, color = "black"),
      legend.position = "top",
      legend.title    = element_blank()
    )
}

# -------------------------------------------------------------------------
# 1. Processing Function
# -------------------------------------------------------------------------
process_bed <- function(file_path, sample_name) {
  read_tsv(file_path, col_names = FALSE, show_col_types = FALSE) |>
    select(
      start   = X2,
      end     = X3,
      ee_id   = X4,
      gene_id = X6
    ) |>
    distinct(gene_id, ee_id, start, end) |>
    count(gene_id, name = "ee_count") |>
    count(ee_count, name = "gene_count") |>
    mutate(sample = sample_name)
}

# -------------------------------------------------------------------------
# 2. Data 
# -------------------------------------------------------------------------
files_to_process <- list(
  list(path = "tomato_EE_list/tomato_permissive_EEs.txt", sample = "Permissive"),
  list(path = "tomato_EE_list/tomato_stringent_EEs.txt",  sample = "Stringent")
)

plot_data <- files_to_process |>
  map_dfr(\(x) process_bed(x$path, x$sample)) |>
  complete(
    ee_count,
    sample,
    fill = list(gene_count = 0)
  )

# -------------------------------------------------------------------------
# 3. Main Plot
# -------------------------------------------------------------------------
main_plot <- ggplot(plot_data, aes(x = factor(ee_count), y = gene_count, fill = sample)) +
  geom_col(
    position = position_dodge(width = 0.8),
    color    = "white",
    width    = 0.7
  ) +
  geom_text(
    aes(label = gene_count),
    position = position_dodge(width = 0.8),
    vjust    = 0.5,
    hjust    = -0.2,
    size     = 3,
    angle    = 90
  ) +
  scale_fill_manual(values = condition_colors) +
  scale_y_continuous(
    limits = c(0, 700),
    expand = expansion(mult = c(0, 0.1))
  ) +
  theme_ee_distribution() +
  labs(
    x = "Number of EEs",
    y = "Number of genes"
  )

print(main_plot)
ggsave("number_of_ees_comparison.svg", main_plot, width = 2.5, height = 4)

# -------------------------------------------------------------------------
# 4. Zoomed Plot (Closeup)
# -------------------------------------------------------------------------
# Filter groups where at least one condition has gene counts below 15
filtered_ee_counts <- plot_data |>
  group_by(ee_count) |>
  filter(any(gene_count < 15)) |>
  pull(ee_count) |>
  unique()

closeup_data <- plot_data |>
  filter(ee_count %in% filtered_ee_counts)

closeup_plot <- ggplot(closeup_data, aes(x = factor(ee_count), y = gene_count, fill = sample)) +
  geom_col(
    position = position_dodge(width = 0.8),
    color    = "white",
    width    = 0.7
  ) +
  geom_text(
    aes(label = gene_count),
    position = position_dodge(width = 0.8),
    vjust    = -0.5,
    size     = 3
  ) +
  scale_fill_manual(values = condition_colors) +
  scale_y_continuous(
    limits = c(0, 150),
    expand = expansion(mult = c(0, 0.1))
  ) +
  theme_ee_distribution() +
  labs(
    x = "Number of EEs",
    y = "Number of genes"
  )

print(closeup_plot)
ggsave("number_of_ees_closeup_comparison.svg", closeup_plot, width = 2.5, height = 4)