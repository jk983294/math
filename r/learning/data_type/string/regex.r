# match
grep("a+", c("abc", "def", "cba a", "aa"), perl = TRUE, value = FALSE)  # [1 3 4] index
grep("a+", c("abc", "def", "cba a", "aa"), perl = TRUE, value = TRUE)  # ['abc', 'cba a', 'aa']
grepl("a+", c("abc", "def", "cba a", "aa"), perl = TRUE, )  # [TRUE, FALSE, TRUE, TRUE]

# find
regexpr("a+", c("abc", "def", "cba a", "aa"), perl = TRUE)  # index [1 -1 3 1] , length [1 -1 1 2]
gregexpr("a+", c("abc", "def", "cba a", "aa"), perl = TRUE)

# extract
x <- c("abc", "def", "cba a", "aa")
regmatches(x, regexpr("a+", x, perl = TRUE))  # ['a', 'a', 'aa']
regmatches(x, gregexpr("a+", x, perl = TRUE))

# replace
sub("(a+)", "z", "abca", perl = TRUE)  # 'zbca' replace first
gsub("(a+)", "z", "abca", perl = TRUE)  # 'zbcz' replace all
# replace reference group with \1~\9
sub("(a+)", "z\\1z", "abca", perl = TRUE)  # 'zazbca' replace first
gsub("(a+)", "z\\1z", "abca", perl = TRUE)  # 'zazbczaz' replace all
# replace with groups
x <- c("abc", "def", "cba a", "aa")
regmatches(x, gregexpr("a+", x, perl = TRUE)) <- list(c("one"), character(0), c("two", 
    "three"), c("four"))

# greedy vs non greedy
number <- "101000000000100"
regmatches(number, gregexpr(pattern = "1.*1", text = number))  # '1010000000001'
regmatches(number, gregexpr(pattern = "1.?1", text = number))  # '101'
