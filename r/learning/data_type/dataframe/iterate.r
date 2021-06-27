dat <- data.frame(x = c(1:5), y = c(11:15))

# by row
apply(dat, 1, function(row) {
  print(paste(row[1], row[2], sep = ","))
})