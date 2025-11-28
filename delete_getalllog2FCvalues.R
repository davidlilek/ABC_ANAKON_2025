# define groups
part1 <- c("C_vs_E", "C_vs_RE", "C_vs_R", "E_vs_R", "E_vs_RE", "R_vs_RE")

# path and names of dataset
ds <- list(
  "L-540 cryo" = "/proj/proteomics/7_combinedHL/evaluation_paper/abc_2025/raw_log2FCdata/hl540_cryo/",
  "L-540 liv"  = "/proj/proteomics/7_combinedHL/evaluation_paper/abc_2025/raw_log2FCdata/hl540_le/",
  "L-428"      = "/proj/proteomics/7_combinedHL/evaluation_paper/abc_2025/raw_log2FCdata/hl428_all/"
)


read_one_set <- function(label, path){
  files <- list.files(path, pattern="\\.csv$", full.names=TRUE)
  if (length(files) == 0) return(NULL)
  dl <- lapply(files, read.csv)
  names(dl) <- part1[seq_along(dl)]
  bind_rows(dl, .id="treatment") |>
    mutate(dataset = label,
           treatment = factor(treatment, levels = part1))
}

final_df <- map2(names(ds), ds, read_one_set) |> list_rbind()

filtered_df <- final_df %>%
  dplyr::filter(adj.P.Val < 0.10,
                abs(logFC) > 2)


write.csv2(filtered_df,"filtered.csv")
