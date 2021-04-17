#!/usr/bin/env Rscript

file.remove("_main.Rmd")

styler::style_dir(".",
  recursive = FALSE,
  filetype = c("R", "Rmd")
)

styler::style_dir("templates",
  recursive = FALSE,
  filetype = c("R", "Rmd")
)

bookdown::render_book("index.Rmd")
