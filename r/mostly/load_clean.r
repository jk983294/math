library(datasets) # Load built-in datasets

summary(iris) # Summary statistics for iris data

search() # currently loaded packages
lsf.str("package:QA") # all exported functions
ls("package:QA") # all exported functions
help(package = QA)

# CLEAN UP #################################################

# Clear packages
detach("package:datasets", unload = TRUE) # For base

# Clear plots
dev.off() # But only if there IS a plot

# Clear console
cat("\014") # ctrl + L