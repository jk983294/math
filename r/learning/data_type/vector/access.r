# Accessing Vector Elements, Indexing starts with position 1
t <- c("Sun", "Mon", "Tue", "Wed", "Thurs", "Fri", "Sat")
t[c(2, 3, 6)]  # ['Mon' 'Tue' 'Fri']
t[c(TRUE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE)]  # ['Sun' 'Fri']
t[c(0, 0, 0, 0, 0, 0, 1)]  # only get index 1 ['Sun']
t[c(-2, -5)]  # drop index 2 and index 5, ['Sun' 'Tue' 'Wed' 'Fri' 'Sat']
t[1]  # ['Sun']
