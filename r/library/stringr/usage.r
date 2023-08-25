library(stringr)

# https://github.com/rstudio/cheatsheets/blob/main/strings.pdf

x <- c("why", "video", "cross", "extra", "deal", "authority")
y <- c("a", "b", "c", "d", "e", "f")

# manage length
str_length(x)
str_pad(x, width = 10L, side = "left", pad = " ")
str_trunc(x, width = 10L, side = "left", pad = " ")
str_trim(c(" x ", " y "), side = "both")
str_squish(c(" x   z ", " y z ")) # trim both ends and collapse multiple spaces into single space

# mutate
str_to_lower(x)
str_to_upper(x)
str_to_title(c("i am good!")) # "I Am Good!"
str_to_sentence(c("i am good!")) # "I am good!"
str_sub(x, start = 1L, end = 2L)
str_sub(x, start = 2L, end = -1L)
str_replace(x, pattern = "[a-e]", replacement = "-") # replace first
str_replace_all(x, pattern = "[a-e]", replacement = "-") # replace all
str_remove(x, pattern = "[a-e]") # remove first
str_remove_all(x, pattern = "[a-e]")

# pattern match
str_detect(x, pattern = "[aeiou]") # return bool list whether it matches pattern
str_starts(x, pattern = "[a-e]")
str_ends(x, pattern = "[a-e]")
str_which(x, pattern = "[a-b]")  # return the indexes which contain the pattern
str_locate(x, pattern = "[aeiou]") # gives the position of the match
str_count(x, pattern = "[aeiou]") # count number of matches in a string

# subset string
str_subset(x, pattern = "[aeiou]") # return only the words match pattern
str_extract(x, pattern = "[aeiou]") # extracts first match
str_extract_all(x, pattern = "[aeiou]") # extracts every match
str_match(x, pattern = "(.)[aeiou](.)") # extracts parts of the match defined by ()
str_match_all(x, pattern = "(.)[aeiou](.)") # extracts parts of the match defined by ()