file.names <- list.files(path = "data/", pattern = "*.csv", full.names = TRUE)
file.list <- lapply(file.names, function(x) {
    read.csv(file = x, header = TRUE)
})

df <- Reduce(function(x, y) {
    merge(x, y, all = TRUE)
}, file.list)
