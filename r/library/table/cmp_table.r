library(data.table)

dt1 <- data.table(
    ID = 1:4,
    B = letters[1:4],
    C = seq(.25, 1, length.out = 4)
)
dt2 <- data.table(
    ID = 4:5,
    B = letters[4:5],
    C = seq(1, 1.25, length.out = 2)
)

# How many of the IDs of dt2 are in dt1
table(dt1$ID %in% dt2$ID)

# Position of IDs of dt1 which appear in dt2
which(dt1$ID %in% dt2$ID)

# IDs of DT1 which appear in dt2
dt1$ID[dt1$ID %in% dt2$ID]

# merge, Get Rid of Duplicate Rows
dt_12 <- merge.data.table(dt1, dt2, all = TRUE)

# Get Number of Exact Duplicate Rows
dt_12b <- rbindlist(list(dt1, dt2)) # row-wise combine
sum(duplicated(dt_12b))

# set operation
union(dt1$ID, dt2$ID)
setdiff(dt1$ID, dt2$ID)
intersect(dt1$ID, dt2$ID)
setequal(dt1$ID, dt2$ID)
is.element(dt1$ID, dt2$ID) # is.element(x, y) identical to x %in% y