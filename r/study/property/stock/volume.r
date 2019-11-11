source("utils/utils.r")
library(ggplot2)

stock_volume_dist_regard_abs_return <- function(symbol) {
    dat <- read_stock_data(symbol)
    df <- dat[, c("close", "volume")]
    names(df) <- c("price", "volume")
    df$return <- diff(log(df[, 1]))
    df <- df[-1, ]  # remove first row
    # cut observation into 3 group based on abs(return)
    df$cuts <- cut(abs(df$return), breaks = c(0, 0.02, 0.04, 0.25), include.lowest = TRUE)
    df$means <- NA  # mean of each group
    for (i in 1:3) {
        group <- which(df$cuts == i)
        if (length(group) > 0) {
            df$means[group] <- mean(df$volume[group])  # calculate each bin's volume mean
        }
    }
    # draw volume dist of return bin
    ggplot(df) + geom_histogram(aes(x = volume)) + facet_grid(cuts ~ .) + geom_vline(aes(xintercept = means), linetype = "dashed", 
        size = 1)
}

stock_volume_dist_regard_abs_return("../data/000001.csv")
