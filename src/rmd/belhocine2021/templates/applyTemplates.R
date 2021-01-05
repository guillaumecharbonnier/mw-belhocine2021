applyTemplateOrlando2014 <- function(orlando,
                                     chunk_label_prefix) {
  src <- knitr::knit_expand("templates/Orlando2014.Rmd")
  res <- knitr::knit_child(
    text = unlist(src),
    envir = environment(),
    quiet = TRUE
  )
  cat(res, sep = "\n")
}

applyTemplateOurChipSeq <- function(gtftk_coverage,
                                    chunk_label_prefix) {
  src <- knitr::knit_expand("templates/OurChipSeq.Rmd")
  res <- knitr::knit_child(
    text = unlist(src),
    envir = environment(),
    quiet = TRUE
  )
  cat(res, sep = "\n")
}
