# elements can be given names and they can be accessed using these names
list_data <- list(c("Jan", "Feb", "Mar"), matrix(c(3, 9, 5, 1, -2, 8), nrow = 2), list("green", 12.3))
names(list_data) <- c("1st Quarter", "A_Matrix", "A Inner list")
print(list_data)

list_data[2]  # access by index
list_data$A_Matrix  # access by name

## single double brackets
my_list <- list(a = c(1, 2), b = matrix(1:10, nrow = 2, ncol = 5), c = data.frame(price = c(89.3), stock = c("MOT")))
my_list[[1]]  # double brackets [[]] return list elements
class(my_list[[1]])  # 'numeric'
my_list[1]  # single brackets return lists
class(my_list[1])  # 'list'
