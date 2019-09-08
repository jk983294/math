as.Date(c("2007-06-22", "2004-02-13"))
as.Date(c("01/05/1965", "08/16/1975"), "%m/%d/%Y")
Sys.Date()  # 返回当天的日期
date()  # 返回当前的日期和时间

# format
today <- Sys.Date()
format(today, format = "%m/%d/%Y")
format(today, format = "%Y-%m-%d")

# Calculations with with dates
startdate <- as.Date("2004-02-13")
enddate <- as.Date("2009-06-22")
enddate - startdate  # Time difference of 1956 days

# Date functions and formatted printing
today <- Sys.Date()
dob <- as.Date("1956-10-12")
difftime(today, dob, units = "weeks")
