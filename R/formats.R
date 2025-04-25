
# %%
read_gmt <- function(file) {
  gmt <- read.table(file, sep = "\t", header = FALSE, fill = NA)
  gs_names <- gmt[, 1]
  genes <- gmt[, -c(1, 2)]
  gs <- split(genes, seq(nrow(genes)))
  gs <- lapply(gs, function(x) unname(t(x))[, 1])
  gs <- lapply(gs, function(x) x[which(x != "")])
  names(gs) <- gs_names
  gs
}

# %%
write_gmt <- function(genesets, fpath) {
  fp <- file(fpath, "w")
  for (n in names(genesets)) {
    set <- genesets[[n]]
    write(glue::glue("{n}\t\t{paste(set, collapse='\t')}"), fp)
  }
  close(fp)
}
