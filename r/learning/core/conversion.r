# vector, array, matrix have same type. list can have different type.  dataframe
# has same type by column

# string > numeric > bool
c("a", TRUE)  # ['a', 'TRUE']
c(3, TRUE)  # [3, 1]

# check
is.numeric(1)
is.character("1")
is.vector(1)
is.matrix(1)
is.data.frame(1)
is.factor(1)
is.logical(1)

# explicit convert
as.numeric("1")
as.character(1)
as.vector(1)
as.matrix(1)
as.data.frame(1)
as.factor(1)
as.logical(1)  # TRUE
