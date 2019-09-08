data <- read.csv("data/input.csv")

# colClasses to avoid string column as factor
data <- read.table("data/input.csv", header = TRUE, row.names = "name", sep = ",", colClasses = c("numeric", "character", "numeric", 
    "character", "character"))
