library(magrittr)

# tee operator %T>%
rnorm(200) %>% matrix(ncol = 2) %T>% plot %>% colSums

# plot() typically don't return anything, use tee operator to skip plot
