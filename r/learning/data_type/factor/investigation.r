data <- c("East", "West", "East", "North", "North", "East", "West", "West", "West", 
    "East", "North")
factor_data <- factor(data)

is.factor(factor_data)  # TRUE
typeof(factor_data)  # 'integer'
attributes(factor_data)  # $levels [1] 'East'  'North' 'West'; $class [1] 'factor'
unclass(factor_data)  # how it stored
as.character(factor_data)  # convert to string vector

# R treats the text column as categorical data and creates factors on it.
height <- c(132, 151, 162, 139, 166, 147, 122)
weight <- c(48, 49, 66, 53, 67, 52, 40)
gender <- c("male", "male", "female", "female", "male", "female", "male")
input_data <- data.frame(height, weight, gender)
is.factor(input_data$gender)
