# remove all variable
rm(list = ls())

# clear packages
detach("package:datasets", unload = TRUE)

# clear console
cat("\014") # ctrl+L