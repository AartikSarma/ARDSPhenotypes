MAST.files <- dir(path = "./output/scMAST/")

for(file in MAST.files) {
  cell <- str_remove(file, " - M.+")
  load(paste("./output/scMAST/", file, sep = ""))
  print(head(fcHurdletest))
  fcHurdletest %>%
    dplyr::select(hgnc = primerid, logFC = coef, fdr) %>%
    filter(fdr < 0.1) %>%
    write_csv(file = paste("./output/exporttoipa/scSeq -", cell, "- MAST.csv"))
}
