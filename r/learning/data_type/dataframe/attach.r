attach(mtcars)  # add dataset mtcars into search path
summary(mpg)  # = summary(mtcars$mpg)
plot(mpg, disp)  # = plot(mtcars$mpg, mtcars$disp)
detach(mtcars)  # remove it from search path

# !!!WARN, the same variable will confilt, the origin one has high priority

# use with can solve same name variable problem
with(mtcars, {
    print(summary(mpg))
    plot(mpg, disp)
    inner_var <- 'x'
    out_var <<- 'y'
})
print(inner_var) # object 'inner_var' not found
print(out_var) # "y"
