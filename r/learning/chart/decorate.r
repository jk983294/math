dose <- c(20, 30, 40, 45, 60)
drugA <- c(16, 20, 27, 40, 60)
drugB <- c(15, 18, 25, 31, 40)

# Adding text, lines, and symbols
plot(dose, drugA, type = "b", col = "red", lty = 2, pch = 2, lwd = 2, main = "Clinical Trials for Drug A", 
    sub = "This is hypothetical data", xlab = "Dosage", ylab = "Drug Response", xlim = c(0, 
        60), ylim = c(0, 70))
lines(dose, drugB, type = "b", pch = 17, lty = 2, col = "blue")
abline(h = c(30), lwd = 1.5, lty = 2, col = "gray")  # 参考线

library(Hmisc)
minor.tick(nx = 3, ny = 3, tick.ratio = 0.5)  # 添加次要刻度线

legend("topleft", inset = 0.05, title = "Drug Type", c("A", "B"), lty = c(1, 2), 
    pch = c(15, 17), col = c("red", "blue"))  # 图例

text(10, 3, "Example of default text")
