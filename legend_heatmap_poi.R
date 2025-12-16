library(ggplot2)

# Farb-Vektor
Category <- c(
  "DNA damage / chromatin / NE stress"          = "#8DD3C7",
  "Mitochondrial dysfunction & energy collapse" = "#FFFFB3",
  "Protein misfolding & ER stress"              = "#BEBADA",
  "Metabolic & stress signaling"                = "#FB8072",
  "Apoptosis"                                   = "#80B1D3",
  "Immune modulation & surface remodeling"      = "#FDB462",
  "Other / assigned to more than one group"     = "#CCCCCC"
)

# Dummy-Plot nur für die Legende
legend_plot <- ggplot(
  data.frame(x = 1:7, y = 1, Category = names(Category)),
  aes(x, y, fill = Category)
) +
  geom_col() +
  scale_fill_manual(
    name   = "Category",   # 👈 Legendentitel
    values = Category
  ) +
  guides(
    fill = guide_legend(nrow = 2)  # 👈 2 Zeilen
  ) +
  theme_void() +
  theme(
    legend.position = "bottom"
  )

legend_plot
