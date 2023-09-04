lines <- readLines("~/github/barn/train/housing.csv")
length(lines)
g <- grep("BAY$", lines)
length(g)
gl <- grepl("BAY$", lines)
length(gl)

# for regex, use stringr package
