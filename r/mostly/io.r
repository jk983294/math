library(data.table)
library(arrow)
library(fst)

in_path_ <- "/home/aaron/junk/y_eval/output/result.feather"
out_path_ <- "/tmp/20210914.1.feather"
dt <- read_feather(in_path_)
write_feather(dt, out_path_)

in_csv_path_ <- "/tmp/20210914.csv"
out_path_ <- "/tmp/20210914.1.fst"
dt_csv <- as.data.table(fread(in_csv_path_))
write.fst(dt_csv, out_path_)

in_path_ <- "/tmp/20210914.fst"
out_path_ <- "/tmp/20210914.1.fst"
dt <- as.data.table(read.fst(in_path_))
write.fst(dt, out_path_)

lines <- readLines("~/github/barn/train/housing.csv")
