library(stringr)

# counting number of characters in a string
nchar("hello")  # 5
str_length("hello")  # 5

# substring
substr("extract", 1, 2)  # 'ex'
x = "extract"
substr(x, 1, 2) <- "EX"  # 'EXtract', replace substr content

# split strings
strsplit(x = "ID-101", split = "-")  #  ['ID', '101']
lets <- strsplit(x = c("ID-101", "ID-102"), split = "-")  # [['ID', '101'], ['ID', '102']]
sapply(lets, "[", 1) # extract first column
sapply(lets, "[", 2) # extract second column
str_split(string = c("ID-101", "ID-102"), pattern = "-", simplify = T)

# find and replace first match
sub(pattern = "L", replacement = "B", x = "ll", ignore.case = T)  # 'Bl'

# find and replace all matches
gsub(pattern = "L", replacement = "B", x = "ll", ignore.case = T)  # 'BB'
gsub(pattern = "\\s", replacement = "", x = "a b  c")  # remove space

# #replace strings
chartr("and", "for", x = " and ")  #letters a,n,d get replaced by f,o,r
str_replace_all(string = "aCityb", pattern = c("City"), replacement = "state")

# get difference between two vectors
setdiff(c("a", "b", "c"), c("a", "c", "f"))  # ['b']

# check if strings are equal
setequal(c("a", "b", "c"), c("a", "c", "f"))  # FALSE
setequal(c("a", "b", "c"), c("a", "c", "b"))  # TRUE

# abbreviate strings
abbreviate(c("monday", "tuesday", "wednesday"), minlength = 3)  # ['mnd', 'tsd', 'wdn']

# upper case, lower case
toupper("hello")  # 'HELLO'
str_to_upper("hello")  # 'HELLO'
tolower("HELLO")  # 'hello'
str_to_lower("HELLO")  # 'hello'
