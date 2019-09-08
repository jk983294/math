# frequency = 12 pegs the data points for every month of a year 
# frequency = 4 pegs the data points for every quarter of a year 
# frequency = 6 pegs the data points for every 10 minutes of an hour 
# frequency = 24*6 pegs the data points for every 10 minutes of a day


# Time Series Analysis
rainfall <- c(799, 1174.8, 865.1, 1334.6, 635.4, 918.5, 685.5, 998.6, 784.2, 985, 
    882.8, 1071)
rainfall.timeseries <- ts(rainfall, start = c(2012, 1), frequency = 12)
print(rainfall.timeseries)
plot(rainfall.timeseries)

# Multiple Time Series
rainfall1 <- c(799, 1174.8, 865.1, 1334.6, 635.4, 918.5, 685.5, 998.6, 784.2, 985, 
    882.8, 1071)
rainfall2 <- c(655, 1306.9, 1323.4, 1172.2, 562.2, 824, 822.4, 1265.5, 799.6, 1105.6, 
    1106.7, 1337.8)
combined.rainfall <- matrix(c(rainfall1, rainfall2), nrow = 12)
rainfall.timeseries <- ts(combined.rainfall, start = c(2012, 1), frequency = 12)
print(rainfall.timeseries)
plot(rainfall.timeseries, main = "Multiple Time Series")
