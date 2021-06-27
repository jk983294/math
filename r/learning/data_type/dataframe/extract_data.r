df <- data.frame(emp_id = c(1:5), emp_name = c("Rick", "Dan", "Michelle", "Ryan", 
    "Gary"), salary = c(623.3, 515.2, 611, 729, 843.25), start_date = as.Date(c("2012-01-01", 
    "2013-09-23", "2014-11-15", "2014-05-11", "2015-03-27")), stringsAsFactors = FALSE)

# Extract Data from data frame
df[1, 2]  # by row and column
df[1, ]  # by row
df[, 2]  # by column
df[["emp_name"]]  # by column
df$emp_name  # by column

all_names <- apply(df, 1, function(x) {
  x[2]  # pure str array, df$emp_name is not pure str vec
})

df[c(TRUE, TRUE, FALSE, FALSE, FALSE), ]  # select by rows
df[df$emp_id < 3, ]  # select by rows with predicate
subset(df, emp_id < 3)  # same with subset
df[df$emp_id < 3 & df$salary > 600, ]  # multi-condition
df[df$emp_name %in% c("Rick", "Dan"), ]  # select by rows with predicate
result <- data.frame(df$emp_name, df$salary)  # extract by column name
result <- df[1:2, ]  # extract first two rows
result <- df[c(3, 5), c(2, 4)]  # extract row (3, 5) and column (2, 4)

# drop some columns
to_drop <- names(df) %in% c("emp_name", "start_date")
df[!to_drop]
