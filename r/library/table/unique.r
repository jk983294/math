library(data.table)
library(dplyr)

dt <- data.table(
    x = c(1, 1, 2, 2, 2, 3, 3, 4, 5),
    group = c("A", "A", "A", "B", "B", "C", "C", "C", "C"),
    group1 = c("a", "a", "b", "b", "b", "b", "c", "c", "c")
)

# base r approach
(ans <- aggregate(dt$x, list(dt$group), function(x) length(unique(x))))
(ans <- aggregate(dt$x, list(dt$group, dt$group1), function(x) length(unique(x))))

# dplyr approach
(ans <- dt %>% group_by(group) %>% dplyr::summarise(count = n_distinct(x)))
(ans <- dt %>% group_by(group, group1) %>% dplyr::summarise(count = n_distinct(x)))

# dt approach
(ans <- dt[ , .(count = length(unique(x))), by = group])
(ans <- dt[ , .(count = length(unique(x))), by = .(group, group1)])
