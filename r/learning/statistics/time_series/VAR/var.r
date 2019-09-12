library(vars)

a <- c(3.1, 3.4, 5.6, 7.5, 5, 6, 6, 7.5, 4.5, 6.7, 9, 10, 8.9, 7.2, 7.5, 6.8, 9.1)
b <- c(3.2, 3.3, 4.6, 8.5, 6, 6.2, 5.9, 7.3, 4, 7.2, 3, 12, 9, 6.5, 5, 6, 7.5)
d <- c(4, 6.2, 5.3, 3.6, 7.5, 6, 6.2, 6.9, 8.2, 4, 4.2, 5, 11, 3, 2.5, 3, 6)
df <- data.frame(a, b, d)

abvar <- VAR(df, p = 2)  #run VAR(2)
coef(abvar)  # coefficient estiamte of VAR
summary(abvar)  # summary of VAR

library(coefplot)

coefplot(abvar$varresult$a)
coefplot(abvar$varresult$b)
coefplot(abvar$varresult$d)
