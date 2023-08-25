name <- "Fred"
age <- 50
anniversary <- as.Date("1991-10-12")
a <- str_glue(
  "My name is {name}, ",
  "my age next year is {age + 1}, ",
  "and my anniversary is {format(anniversary, '%A, %B %d, %Y')}."
)

str_glue("My name is {name}, not {{name}}.")

# interpolating string
mtcars %>% str_glue_data("{rownames(.)} has {hp} hp")
