library(tidyverse)

raw = read_lines("~/Desktop/ast.txt") |>
    paste0(collapse=" ") |>
    str_squish() |>
    str_replace_all("[<>]", ".") |>
    str_replace_all(" \\(", "(")

cat(str_sub(raw, 1, 500), "\n")

