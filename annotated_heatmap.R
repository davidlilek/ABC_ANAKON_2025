## -------------------------------
## 1) Datensätze einlesen & vorbereiten
## -------------------------------

# Vergleichsgruppen
part1 <- c("C_vs_E", "C_vs_RE", "C_vs_R", "E_vs_R", "E_vs_RE", "R_vs_RE")

# Pfade zu den Datensätzen
ds <- list(
  "L-540 cryo" = "/proj/proteomics/7_combinedHL/evaluation_paper/abc_2025/raw_log2FCdata/hl540_cryo/",
  "L-540 liv"  = "/proj/proteomics/7_combinedHL/evaluation_paper/abc_2025/raw_log2FCdata/hl540_le/",
  "L-428"      = "/proj/proteomics/7_combinedHL/evaluation_paper/abc_2025/raw_log2FCdata/hl428_all/"
)

# Hilfsfunktion: ein Set (Ordner) einlesen
read_one_set <- function(label, path){
  files <- list.files(path, pattern = "\\.csv$", full.names = TRUE)
  if (length(files) == 0) return(NULL)
  dl <- lapply(files, read.csv)
  names(dl) <- part1[seq_along(dl)]
  d <- dplyr::bind_rows(dl, .id = "treatment")
  
  d |>
    dplyr::filter(Gene.name %in% proteins_of_interest) |>
    dplyr::mutate(
      dataset   = label,
      treatment = factor(treatment, levels = part1)
    )
}

# alle Datensätze einlesen und zusammenführen
final_df <- purrr::map2(names(ds), ds, read_one_set) |> purrr::list_rbind()

# logFC begrenzen
final_df$logFC[!is.na(final_df$logFC) & final_df$logFC >  4] <-  4
final_df$logFC[!is.na(final_df$logFC) & final_df$logFC < -4] <- -4

# Spalten-ID (für Matrix) + Anzeigename (nur Treatment)
final_df <- final_df |>
  dplyr::mutate(
    col_id      = paste0(dataset, " | ", treatment),
    col_display = treatment
  )

## -------------------------------
## 2) logFC-Matrix & Signifikanz-Matrix
## -------------------------------

# logFC-Matrix
mat_fc <- final_df |>
  dplyr::select(Gene.name, col_id, logFC, dataset, treatment, col_display) |>
  dplyr::distinct() |>
  dplyr::arrange(dataset, treatment) |>
  dplyr::select(Gene.name, col_id, logFC) |>
  tidyr::pivot_wider(names_from = col_id, values_from = logFC) |>
  tibble::column_to_rownames("Gene.name") |>
  as.matrix()

mat_fc[is.na(mat_fc)] <- 0

# Signifikanz-Matrix (Sternchen)
mat_sig <- final_df |>
  dplyr::mutate(
    Significant = dplyr::case_when(
      adj.P.Val < 0.001 ~ "\u2731\u2731\u2731",
      adj.P.Val < 0.01  ~ "\u2731\u2731",
      adj.P.Val < 0.05  ~ "\u2731",
      adj.P.Val < 0.1   ~ ".",
      TRUE              ~ " "
    )
  ) |>
  dplyr::select(Gene.name, col_id, Significant, dataset, treatment) |>
  dplyr::distinct() |>
  dplyr::arrange(dataset, treatment) |>
  dplyr::select(Gene.name, col_id, Significant) |>
  tidyr::pivot_wider(names_from = col_id, values_from = Significant) |>
  tibble::column_to_rownames("Gene.name") |>
  as.matrix()

mat_sig[is.na(mat_sig)] <- " "

# Spalten-Annotation (Datasets)
annot_col <- final_df |>
  dplyr::distinct(col_id, dataset, treatment, col_display) |>
  dplyr::arrange(dataset, treatment) |>
  tibble::column_to_rownames("col_id")

## -------------------------------
## 3) Zeilen-Annotation: funktionelle Gruppen
## -------------------------------

# Mapping Gen -> Kategorie (nach deinem Text)
gene2cat <- c(
  # DNA Damage, Chromatin & Nuclear Envelope Stress
  "APEX1" = "DNA damage / chromatin / NE stress",
  "ARMT1" = "DNA damage / chromatin / NE stress",
  "BANF1" = "DNA damage / chromatin / NE stress",
  "BCCIP" = "DNA damage / chromatin / NE stress",
  "PARP1" = "DNA damage / chromatin / NE stress",
  
  # Mitochondrial Dysfunction & Energy Metabolism
  "CYC1"   = "Mitochondria & energy metabolism",
  "SDHA"   = "Mitochondria & energy metabolism",
  "SDHB"   = "Mitochondria & energy metabolism",
  "OGDH"   = "Mitochondria & energy metabolism",
  "DLST"   = "Mitochondria & energy metabolism",
  "SUCLG1" = "Mitochondria & energy metabolism",
  "SOD2"   = "Mitochondria & energy metabolism",
  "VDAC1"  = "Mitochondria & energy metabolism",
  "VDAC2"  = "Mitochondria & energy metabolism",
  "VDAC3"  = "Mitochondria & energy metabolism",
  "NENF"   = "Mitochondria & energy metabolism",
  
  # 3.3.3 Protein Misfolding, ER Stress & Proteostasis
  "BAG6"     = "Protein misfolding & ER stress",
  "SERPINB6" = "Protein misfolding & ER stress",
  "CDC37"    = "Protein misfolding & ER stress",
  "NUDC"     = "Protein misfolding & ER stress",
  "MANF"     = "Protein misfolding & ER stress",
  "MYDGF"    = "Protein misfolding & ER stress",
  "PDCD6"    = "Protein misfolding & ER stress",
  "PDCD6IP"  = "Protein misfolding & ER stress",
  "PCYT1A"   = "Protein misfolding & ER stress",
  
  # 3.3.4 Metabolic & Stress Signaling Regulation
  "HK1"    = "Metabolic & stress signaling",
  "HK2"    = "Metabolic & stress signaling",
  "STAT1"  = "Metabolic & stress signaling",
  "STAT3"  = "Metabolic & stress signaling",
  "STAT5A" = "Metabolic & stress signaling",
  "STAT6"  = "Metabolic & stress signaling",
  "GRB2"   = "Metabolic & stress signaling",
  "PURA"   = "Metabolic & stress signaling",
  "NUP214" = "Metabolic & stress signaling",
  "OGFR"   = "Metabolic & stress signaling",
  "RBM4"   = "Metabolic & stress signaling",
  "NUDT5"  = "Metabolic & stress signaling",
  "SLC4A7" = "Metabolic & stress signaling",
  "KCNAB2" = "Metabolic & stress signaling",
  
  # 3.3.5 Apoptosis
  "AIFM1"   = "Apoptosis",
  "CASP3"   = "Apoptosis",
  "LGALS1"  = "Apoptosis",
  "LGALS3"  = "Apoptosis",
  "PDCD5"   = "Apoptosis",
  "PRCP"    = "Apoptosis",
  "TNFAIP8" = "Apoptosis",
  "TNFRSF8" = "Apoptosis",
  
  # 3.3.6 Immune Modulation & Surface Remodeling
  "CD63" = "Immune modulation & surface remodeling",
  "CD70" = "Immune modulation & surface remodeling",
  "HDGF" = "Immune modulation & surface remodeling",
  "LGALS1" = "Immune modulation & surface remodeling",
  "LGALS3 " = "Immune modulation & surface remodeling"
)

# alle Gene (Zeilen) der Matrix
all_genes <- rownames(mat_fc)

# Kategorie-Vektor; Standard: "Other / not assigned"
row_cat <- gene2cat[all_genes]
row_cat[is.na(row_cat)] <- "Other / not assigned"

annotation_row <- data.frame(
  Category = factor(
    row_cat,
    levels = c(
      "DNA damage / chromatin / NE stress",
      "Mitochondria & energy metabolism",
      "Protein misfolding & ER stress",
      "Metabolic & stress signaling",
      "Apoptosis",
      "Immune modulation & surface remodeling",
      "Other / not assigned"
    )
  ),
  row.names = all_genes
)

# Zeilen nach Kategorie + Genname sortieren
ord <- order(annotation_row$Category, rownames(mat_fc))
mat_fc         <- mat_fc[ord, , drop = FALSE]
mat_sig        <- mat_sig[ord, , drop = FALSE]
annotation_row <- annotation_row[ord, , drop = FALSE]

# gaps_row: dort, wo Kategorie wechselt
cat_vec  <- as.character(annotation_row$Category)
gaps_row <- which(head(cat_vec, -1) != tail(cat_vec, -1))

## -------------------------------
## 4) Farben & pheatmap-Aufruf
## -------------------------------

# Farbskala für logFC
rng    <- 4
breaks <- seq(-rng, rng, length.out = 101)
col_palette <- colorRampPalette(c("blue", "white", "red"))(length(breaks))

# Spalten-Gruppierung (Datasets)
col_gaps <- c(6, 12)

# Annotation-Farben
anno_cols <- list(
  dataset = c(
    "L-540 cryo" = "#D9D9D9",
    "L-540 liv"  = "#D9D9D9",
    "L-428"      = "#D9D9D9"
  ),
  Category = c(
    "DNA damage / chromatin / NE stress"             = "#8DD3C7",
    "Mitochondria & energy metabolism"               = "#FFFFB3",
    "Protein misfolding & ER stress"                 = "#BEBADA",
    "Metabolic & stress signaling"                   = "#FB8072",
    "Apoptosis"                                      = "#80B1D3",
    "Immune modulation & surface remodeling"         = "#FDB462",
    "Other / not assigned"                           = "#CCCCCC"
  )
)

# Heatmap
pheatmap::pheatmap(
  mat_fc,
  color = col_palette,
  breaks = breaks,
  border_color = NA,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = mat_sig,
  number_color = "#00000080",
  fontsize_number = 3,
  show_rownames = TRUE,
  show_colnames = TRUE,
  labels_col = annot_col$col_display,
  annotation_col  = annot_col[, c("dataset"), drop = FALSE],
  annotation_row  = annotation_row,
  annotation_colors = anno_cols,
  annotation_legend = TRUE,
  annotation_names_col = TRUE,
  gaps_col = col_gaps,
  gaps_row = gaps_row,
  fontsize_row = 4.65,
  fontsize_col = 6,
  cellheight = 3.75,
  cellwidth  = 12,
  fontsize   = 6
)
