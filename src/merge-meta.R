#!/usr/bin/env Rscript

library(optparse)
library(digest)

option_list <- list(
    make_option(c("-c", "--clusters"), default="",
                help="Clustering output"),
    make_option(c("-m", "--metadata"), default="",
                help="Metadata including PLZ"),
    make_option(c("-o", "--outfile"), default="",
                help="Output file (tsv) with clustering and metadata merged"),
    make_option(c("-d", "--seed"), default=531,
                help="Set the random seed at the beginning of the script")
    )
args <- parse_args(OptionParser(option_list = option_list), convert_hyphens_to_underscores = TRUE)
conf <- args
conf$id <- substr(digest(args), 1, 8)
set.seed(args$seed)

library(dplyr, warn.conflicts = FALSE)
library(readr, warn.conflicts = FALSE)

clust <- read.delim(args$clusters)
meta <- read.delim(args$metadata)

metaclust <- clust %>% left_join(meta, by = c("seqName" = "ID"))

write_tsv(metaclust, args$outfile)

