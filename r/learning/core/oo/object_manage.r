# find variables
ls()
ls(pattern = "var")  # regex find variable name
ls(all.name = TRUE)  # variables starting with dot(.) are hidden

# deleting Variables
var.3 = 5
rm(var.3)
print(var.3)  # object 'var.3' not found
rm(list = ls())  # rm all variables
print(ls())  # character(0)

# save
z <- rnorm(1000)
save(z, file="zfile")
rm(z)
exists("z") # FALSE

# load
load(file="zfile")
exists("z") # TRUE
