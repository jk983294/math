# elements can be given names and they can be accessed using these names
list_data <- list(c("Jan", "Feb", "Mar"), matrix(c(3, 9, 5, 1, -2, 8), nrow = 2), 
    list("green", 12.3))
names(list_data) <- c("1st Quarter", "A_Matrix", "A Inner list")
print(list_data)

list_data[2]  # access by index
list_data$A_Matrix  # access by name
