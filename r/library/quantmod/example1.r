library(quantmod)

getSymbols("AAPL", src = "yahoo")
head(AAPL)
barChart(AAPL)
# Add multi-coloring and change background to white
candleChart(AAPL, multi.col = TRUE, theme = "white")
