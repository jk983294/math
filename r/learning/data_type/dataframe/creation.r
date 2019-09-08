# data frame, tabular data objects, each column can contain different modes of
# data, a list of vectors of equal length.
bmi <- data.frame(gender = c("Male", "Male", "Female"), height = c(152, 171.5, 165), 
    weight = c(81, 93, 78), Age = c(42, 38, 26))
colnames(bmi) <- c("Gender", "Height", "Weight", "Age")
names(bmi)[2] <- "height"  # change column name

# concat several columns into a data frame
emp.data <- data.frame(emp_id = c(1:5), emp_name = c("Rick", "Dan", "Michelle", "Ryan", 
    "Gary"), salary = c(623.3, 515.2, 611, 729, 843.25), start_date = as.Date(c("2012-01-01", 
    "2013-09-23", "2014-11-15", "2014-05-11", "2015-03-27")), stringsAsFactors = FALSE)
