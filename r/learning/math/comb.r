require(combinat)

get_all_comb <- function(seq) {
    b <- list()
    for (i in seq) {
        b <- c(b, list(0:(i - 1)))
    }
    expand.grid(b)
}

# all combination with repetitions
expand.grid(list(1:2, 1:3)) # 第一个元素2个选择,第二个元素3个选择
get_all_comb(c(2, 3))
permn(3) # A(3, 3)
permn(0:2) # 每个元素3个选择

# all combination without repetitions
combn(3, 2)
combn(1:6, 2) # 6 choose 2
choose(6, 2) # only give number, not combinations