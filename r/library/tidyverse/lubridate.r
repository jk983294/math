library(lubridate)

# https://rawgit.com/rstudio/cheatsheets/main/lubridate.pdf

# parse
d1 <- lubridate::ymd(20101215L)
lubridate::ymd_hms("2017-11-28T14:02:00")
lubridate::ymd_hms("2017-11-28 14:02:00")
lubridate::mdy("4/1/17")
lubridate::yq("2001:Q3")

# query component
lubridate::year(d1)
lubridate::month(d1)
lubridate::day(d1)
lubridate::hour(d1)
lubridate::minute(d1)
lubridate::second(d1)
lubridate::wday(d1)
lubridate::wday(d1, label = TRUE)

# period
p1 <- lubridate::period(5L, unit = "days")
lubridate::period_to_seconds(p1)

d2 <- lubridate::ymd(20101217L)
p2 <- lubridate::as.period(d2 - d1)
lubridate::period_to_seconds(p2)

# duration
(x1 <- lubridate::duration(5L, unit = "days"))
