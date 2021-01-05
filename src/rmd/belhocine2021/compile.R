#!/usr/bin/env Rscript

styler::style_dir(".",
  recursive = FALSE,
  filetype = c("R", "Rmd")
)

styler::style_dir("templates",
  recursive = FALSE,
  filetype = c("R", "Rmd")
)

file.remove("_main.Rmd")
bookdown::render_book("index.Rmd")
