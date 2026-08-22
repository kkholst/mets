# Rebuild pre-built vignettes from .Rmd.orig sources
# Run from the vignettes/ directory
#
# Usage:
#   Rscript rebuild-vignettes.R              # rebuild all
#   Rscript rebuild-vignettes.R foo.Rmd.orig # rebuild one

args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 0) {
  files <- args
} else {
  files <- list.files(".", pattern = "\\.Rmd\\.orig$", full.names = TRUE)
}

for (f in files) {
  outfile <- sub("\\.orig$", "", f)
  vignette_name <- tools::file_path_sans_ext(tools::file_path_sans_ext(basename(f)))
  message("Knitting ", f, " -> ", outfile)
  knitr::opts_chunk$set(
    dev = "png",
    dev.args = list(type = "cairo"),
    dpi = 96,
    fig.cap = "",
    fig.path = paste0("figure/", vignette_name, "-")
  )
  knitr::knit(f, output = outfile)
  # Prepend auto-generated comment after YAML header
  lines <- readLines(outfile)
  yaml_markers <- which(lines == "---")
  if (length(yaml_markers) >= 2) {
    yaml_end <- yaml_markers[2]
    lines <- append(lines,
      "<!-- DO NOT EDIT: auto-generated from the .Rmd.orig source -->",
      after = yaml_end)
  }
  writeLines(lines, outfile)
  message("Done: ", outfile)
}
