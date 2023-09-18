library(datasets) # Load built-in datasets
# later loaded package appear front in search env list
# when search a symbol, first ".GlobalEnv" then user loaded packages, at last package:base
# check it with search()
library(QA)

summary(iris) # Summary statistics for iris data

# below LD_LIBRARY_PATH can avoid some install.packages error
Sys.setenv(LD_LIBRARY_PATH=paste("/opt/anaconda3/lib", Sys.getenv("LD_LIBRARY_PATH"),sep=":"))
install.packages("arrow")
remove.packages("arrow")

update.packages()  # update all available packages
update.packages(ask = FALSE) # update, without prompts

search() # currently loaded packages
lsf.str("package:QA") # all exported functions
ls("package:QA") # all exported functions
help(package = QA)
sessionInfo() # check pkg version
packageVersion("data.table")

# query function
str(lm)
args(lm)
formals(lm)  # argument list
body(lm)

# CLEAN UP #################################################

# Clear packages
detach("package:datasets", unload = TRUE) # For base

# Clear plots
dev.off() # But only if there IS a plot

# Clear console
cat("\014") # ctrl + L