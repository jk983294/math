library(dplyr)
library(data.table)

data <- data.table(x1 = c("a", "b", "c", "b"), x2 = c(2L, 1L, 3L, 4L), x3 = c(1L, 3L, 2L, 2L))

data[order(x2), ]
data[order(x1, -x3)] # multi dimension sort
data[order(-x2), ] # reverse order
dplyr::arrange(data, x2) # order data with dplyr

# key
(setorder(data, x1, -x2))
(setkey(data, x1, x2))
(setkeyv(data, c("x1", "x2")))
key(data) # returns all the key columns

# secondary indices
# Why do we need secondary indices? – Reordering a data.table can be expensive and not always ideal
(setindex(data, x2)) # set secondary index for further calculation
(setindexv(data, c("x2"))) # set secondary index for further calculation
indices(data) # returns all indices

# auto create secondary indices when we use == or %in% on a single column for the first time
(ans <- data[x3 == 2L])
(ans <- data[x3 %in% 1L:2L])

# subset by key
data[.("b", 1L)]
data[.("b")]
data[J("b")] # J : it is a direct alias of list
data[.("b"), mult = "first"] # first match
data[.("b"), mult = "last"] # last match
data[.(c("a", "b"))]
data[.(unique(x1), 1L)]

# subset on indices
data[.(2L), on = "x3"]
data[.(3L, 2L), on = c("x2", "x3")]
data[.(2L), on = "x3", mult = "first"] # first match
data[.(2L), on = "x3", mult = "last"] # last match

# subset calc
data[.("b"), max(x3)]
data[.(2L), max(x2), on = "x3"]

# sub-assign
data[.("b"), x3 := x3 + 1L]
data[.(2L), x2 := x2 + 1L, on = "x3"]

# group
data[.("b"), max(x3), by = .(x2)]
data[.(2L), max(x2), by = .(x2), on = "x3"]
