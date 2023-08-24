library(datasets)

head(iris$Sepal.Length, 5)
tail(iris$Sepal.Length, 5)

summary(iris$Species) # Categorical variable
summary(iris$Sepal.Length) # Quantitative variable
psych::describe(iris$Sepal.Length)

# handle str vec
(vec <- c("A", "B", "A", "C", "C", "A", "B", "C", "D"))
(f_vec <- as.factor(vec))
(dt <- data.table(f_vec))
(dplyr::count(dt, f_vec))  # Get Frequency
dt[ , .N, by = f_vec] # Get Frequency
