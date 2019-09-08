# Need some data to play with
df1 <- data.frame(LETTERS, dfindex = 1:26)
df2 <- data.frame(letters, dfindex = c(1:10, 15, 20, 22:35))

# INNER JOIN: returns rows when there is a match in both tables.
merge(df1, df2)
merge(df1, df2, by="dfindex")

# FULL (outer) JOIN: all records from both the tables and fill in NULLs for
# missing matches on either side.
merge(df1, df2, all = TRUE)


# join by name
names(df1) <- c("alpha", "lotsaNumbers")
merge(df1, df2, by.x = "lotsaNumbers", by.y = "dfindex")

# Outer join: merge(x = df1, y = df2, by = 'CustomerId', all = TRUE)

# Left outer: merge(x = df1, y = df2, by = 'CustomerId', all.x = TRUE)

# Right outer: merge(x = df1, y = df2, by = 'CustomerId', all.y = TRUE)

# Cross join: merge(x = df1, y = df2, by = NULL)

# merge on multiple columns: merge(df1, df2, by = c('CustomerId', 'OrderId'))
