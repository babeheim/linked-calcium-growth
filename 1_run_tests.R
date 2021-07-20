
rm(list = ls())

source("project_support.R")

test_files <- list.files("./tests", pattern = "^test_.*R$",
  full.names = TRUE, recursive = FALSE)

for (i in seq_along(test_files)) source(test_files[i])
