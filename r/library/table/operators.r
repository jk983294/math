# dt[i, j, by]
# take dt, subset or reorder rows using i, calculate j, grouped by by
# i is for operations you’d do on rows (choosing rows based on row numbers or conditions);
# j is what you’d do with columns (select columns or create new columns from calculations).

help("special-symbols")

dt <- data.table(x = c("a", "b", "a", "b"), y = c(2., 1., 3., NA), z = c(1., 3., NA, 2.), m = c("m", "m", "n", "m"))

# walrus := operator, assignment by reference operator
# because it changed the existing object by adding the new column without needing to save it to a new variable.
dt[, has_a := ifelse(x %like% "a", TRUE, FALSE)]
dt[, has_m := ifelse(m %like% "\\bm\\b", TRUE, FALSE)] # \\b means a word boundary

# add multi columns at once, turn that walrus operator into a function by backquoting it
dt[, `:=`(
    has_a = ifelse(x %like% "a", TRUE, FALSE),
    has_m = ifelse(m %like% "\\bm\\b", TRUE, FALSE)
)]

# %between% operator
dt[y %between% c(1.5, 3.5)]

# fcase() function is similar to SQL’s CASE WHEN statement and dplyr’s case_when() function
# The basic syntax is fcase(condition1, "value1", condition2, "value2")
dt[, am_comb := fcase(
    has_a & !has_m, "a",
    !has_a & has_m, "m",
    has_a & has_m, "Both",
    !has_a & !has_m, "Neither"
)]