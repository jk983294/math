bmi <- data.frame(gender = c("Male", "Male", "Female"), height = c(152, 171.5, 165), 
    weight = c(81, 93, 78), Age = c(42, 38, 26))
colnames(bmi) <- c("Gender", "Height", "Weight", "Age")

names(bmi)[2] <- "height"  # change column name
