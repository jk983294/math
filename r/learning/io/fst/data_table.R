library(fst)

# Generate some random data frame with 10 million rows and various column types
nr_of_rows <- 1e+07

df <- data.frame(Logical = sample(c(TRUE, FALSE, NA), prob = c(0.85, 0.1, 0.05), nr_of_rows, replace = TRUE), Integer = sample(1L:100L, 
    nr_of_rows, replace = TRUE), Real = sample(sample(1:10000, 20)/100, nr_of_rows, replace = TRUE), Factor = as.factor(sample(labels(UScitiesD), 
    nr_of_rows, replace = TRUE)))

# Store the data frame to disk
write.fst(df, "/tmp/dataset.fst")
# write.csv(df, '/tmp/dataset.csv') # 压缩比大概是6:1

# Retrieve the data frame again
df1 <- read.fst("/tmp/dataset.fst")
