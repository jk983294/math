library(magrittr)

# unary function
f1 <- . %>% cos() %>% sin()
f1
f1(0)

f2 <- function(.) sin(cos(.))
f2
f2(0)
