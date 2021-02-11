library("data.table")
library("rbgen")
set.seed(1234)

# load in summary stats
gwas <- fread("output.tsv")

# load in BGEN file
ranges <- data.frame(
  chromosome = "01",
  start = 0,
  end = .Machine$integer.max
)
snp <- bgen.load("../../lib/bgen/example/example.v11.bgen", ranges)