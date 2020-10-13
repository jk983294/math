# RStudio usage
# Ctrl + Enter              line by line execute
# Ctrl + L                  clear console

# run script
Rscript hello_world.r
Rscript -e 'install.packages(c("devtools", "Rcpp"))'

# no such file or directory
sudo chown -R $(whoami) ~/.rstudio*

# batch mode
R CMD BATCH options to_execute.R out.Rout
R CMD BATCH hello_world.r /tmp/out.log

# turn off “Hit <Return> to see next plot”
par(ask=F)
devAskNewPage(ask = FALSE)
