first_last <- function(...) {
    arg1 <- ..1
    n <- ...length()
    names <- ...names()
    last_arg <- ...elt(n)
    a <- list(arg1, last_arg)
    names(a) <- c(names[1], names[n])
    a
}

first_last(x = 1, y = 2, z = 3, m = 4)

create_df <- function(...) {
    dots <- list(...)
    dots_names <- names(dots)

    if (is.null(dots_names) || "" %in% dots_names) {
        stop("all arguments must be named")
    }

    list2DF(dots)
}

create_df(1L:10L, letters[1L:10L])
create_df(x = 1L:10L, y = letters[1L:10L])
