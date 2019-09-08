format(23.123456789, digits = 4)  # '23.12'
format(c(6, 13.14521), scientific = TRUE)  # '6.000000e+00' '1.314521e+01'
format(23.47, nsmall = 5)  # '23.47000'
format(13.7, width = 6)  # '  13.7'
format("Hello", width = 8, justify = "l")  # 'Hello   '
format("Hello", width = 8, justify = "c")  # ' Hello  '
format("Hello", width = 8, justify = "r")  # '   Hello'

# sprintf
sprintf("data %d", 1)

# paste
a <- "hello"
b <- "world"
paste(a, b)  # 'hello world'
paste(a, b, sep = "-")  # 'hello-world'
paste(a, b, 1:2, sep = "-", collapse = "*")  # 'hello-world-1*hello-world-2'

# cat
cat("hello", "world", sep = " ")  # hello world
cat("hello", "world", "!")  # hello world !

# nice look
string <- "Los Angeles, officially the City of Los Angeles and often known by its initials L.A., is the second-most populous city in the United States (after New York City), the most populous city in California and the county seat of Los Angeles County. Situated in Southern California, Los Angeles is known for its Mediterranean climate, ethnic diversity, sprawling metropolis, and as a major center of the American entertainment industry."
strwrap(string)
