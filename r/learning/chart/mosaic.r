ftable(Titanic)
library(vcd)
mosaic(Titanic, shade = TRUE, legend = TRUE)

library(vcd)
mosaic(~Class + Sex + Age + Survived, data = Titanic, shade = TRUE, legend = TRUE)
