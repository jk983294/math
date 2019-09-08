emp <- data.frame(emp_id = c(1:5), emp_name = c("Rick", "Dan", "Michelle", "Ryan", 
    "Gary"), salary = c(623.3, 515.2, 611, 729, 843.25), start_date = as.Date(c("2012-01-01", 
    "2013-09-23", "2014-11-15", "2014-05-11", "2015-03-27")), stringsAsFactors = FALSE)

# get the structure of the Data Frame
typeof(emp)  # 'list'
str(emp)

# summary of data in data frame
summary(emp)

dim(emp)  # row column
nrow(emp)  # row
ncol(emp)  # column

# show part data
head(object)  # 列出某个对象的开始部分
tail(object)  # 列出某个对象的最后部分
