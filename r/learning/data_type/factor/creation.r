# Factors are used to categorize the data and store it as levels it can store
# number and string
data <- c("East", "West", "East", "North", "North", "East", "West", "West", "West", 
    "East", "North")
factor_data <- factor(data)  # default order is based on alphabeta

# Changing the order of Levels
new_order_data <- factor(factor_data, levels = c("East", "West", "North"))  # override default order by providing levels

# Generating Factor Levels
v <- gl(2, 2, labels = c("Tampa", "Seattle"))  # [1] Tampa   Tampa   Seattle Seattle,  Levels: Tampa Seattle

# ordered factor
status <- c("Poor", "Improved", "Excellent", "Poor")
status <- factor(status, ordered = TRUE)  # stored as (3, 2, 1, 3), it will auto choose correct analysis based on ordered fact

# create un-ordered factor from number given by level and label
sex <- c(1, 1, 2, 2)
sex <- factor(sex, levels = c(1, 2), labels = c("Male", "Female"))
