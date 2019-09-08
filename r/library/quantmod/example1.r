getSymbols("AAPL", src = "yahoo")
barChart(AAPL)
# Add multi-coloring and change background to white
candleChart(AAPL, multi.col = TRUE, theme = "white")
