# creating a time series object
sales <- c(18, 33, 41, 7, 34, 35, 24, 25, 24, 21, 25, 20, 22, 31, 40, 29, 25, 21, 
    22, 54, 31, 25, 26, 35)
tsales <- ts(sales, start = c(2003, 1), frequency = 12)
tsales
plot(tsales)

start(tsales) # 2003 1
end(tsales) # 2004 12
frequency(tsales) # 12

tsales.subset <- window(tsales, start = c(2003, 5), end = c(2004, 6))
tsales.subset
