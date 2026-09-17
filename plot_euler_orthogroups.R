library(tidyverse)
library(eulerr)

# -------------------------------------------------------------------------
# 0. Global Setup & Color Palette
# -------------------------------------------------------------------------
set_colors <- c(
  "EEs Arabidopsis"    = "#78C8F8FF",
  "Tomato Permissive"  = "#DD88ACFF",
  "Tomato Stringent"   = "#9E4770FF"
)

# -------------------------------------------------------------------------
# 1. Load Data & Prepare Sets
# -------------------------------------------------------------------------
permissive_orthogroups <- read_csv2("EEs_Orthogroups_permissive.csv", show_col_types = FALSE)
stringent_orthogroups  <- read_csv2("EEs_Orthogroups_stringent.csv", show_col_types = FALSE)

euler_set_list <- list(
  "EEs Arabidopsis"   = permissive_orthogroups |> filter(EEs_arabidopsis) |> pull(Orthogroups) |> na.omit(),
  "Tomato Permissive" = permissive_orthogroups |> filter(ATAC_tomato)     |> pull(Orthogroups) |> na.omit(),
  "Tomato Stringent"  = stringent_orthogroups  |> filter(ATAC_tomato)     |> pull(Orthogroups) |> na.omit()
)

# -------------------------------------------------------------------------
# 2. Fit Proportional Euler Model
# -------------------------------------------------------------------------
euler_fit <- euler(euler_set_list, shape = "ellipse")

# -------------------------------------------------------------------------
# 3. Euler Diagram Visualization & Export
# -------------------------------------------------------------------------
euler_plot <- plot(
  euler_fit,
  quantities = list(type = "counts", cex = 1.2, font = 2, col = "gray30"),
  fills      = list(
    fill  = unname(set_colors),
    alpha = 0.25
  ),
  edges      = list(
    col = unname(set_colors),
    lex = 3.5
  ),
  labels     = FALSE,
  legend     = list(side = "top", font = 1, cex = 1)
)

print(euler_plot)
ggsave("euler_arabidopsis_tomato_stringent_permissive.svg", plot = euler_plot, width = 5.5, height = 5.5)