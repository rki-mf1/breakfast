#!/usr/bin/env Rscript

library(optparse)
library(digest)

option_list <- list(
    make_option(c("-i", "--input"), default="",
                help="File with clustering and PLZ, SMAPLING_DATE, etc [default %default]"),
    make_option(c("-o", "--outdir"), default="output/",
                help="Output directory [default %default]")
    )
args <- parse_args(OptionParser(option_list = option_list), convert_hyphens_to_underscores = TRUE)

library(dplyr, warn.conflicts = FALSE)
library(lubridate, warn.conflicts = FALSE)
library(readr, warn.conflicts = FALSE)

if (!dir.exists(args$outdir))
  dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

x <- read.delim(args$input)
#x %>% filter(!is.na(cluster_id)) %>% arrange(cluster_id, desc(SAMPLING_DATE)) %>% head()
cluster_dates <- x %>%
  mutate(SAMPLING_DATE2 = ymd(SAMPLING_DATE)) %>%
  group_by(cluster_id) %>%
  summarize(min = min(SAMPLING_DATE2), max = max(SAMPLING_DATE2), date_range = time_length(interval(min(SAMPLING_DATE2), max(SAMPLING_DATE2)), "day"), n = n(), frac_top_plz = max(table(PRIMARY_DIAGNOSTIC_LAB_PLZ)) / n()) %>%
  mutate(score = (n + 1)/ (date_range + 1)) %>%
  arrange(desc(score))

print(cluster_dates)

# Recent only
recent_cluster_dates <- cluster_dates %>% filter(max >= today() - weeks(2))

write_tsv(cluster_dates, file.path(args$outdir, "cluster-summary.tsv"))
write_tsv(recent_cluster_dates, file.path(args$outdir, "cluster-summary-last-week.tsv"))
