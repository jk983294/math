x <- 5L
if (x > 5) {
    print("x > 5")
} else if (x == 5) {
    print("x is 5")
} else {
    print("x < 5")
}

x <- c("what", "is", "truth")
if ("Truth" %in% x) {
    print("Truth is found")
} else {
    print("Truth is not found")
}

score <- 0.7
outcome <- ifelse(score > 0.5, "Passed", "Failed")
