library(magrittr)
x <- c(0.109, 0.359, 0.63, 0.996, 0.515, 0.142, 0.017, 0.829, 0.907)

x %>% log()  # log(x)
x %>% log() %>% diff() %>% exp() %>% round(1)  # round(exp(diff(log(x))), 1)

# dot placeholder, default is first argument
"Ceci n'est pas une pipe" %>% gsub("une", "un", .)
6 %>% round(pi, digits = .)

# placehold pitfall
ma <- matrix(1:12, 3, 4)
ma %>% max(nrow(ma), ncol(ma))  # 12 max(ma, nrow(ma), ncol(ma))
ma %>% {
    max(nrow(.), ncol(.))
}  # 4
