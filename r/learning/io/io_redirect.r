source("../hello_world.r")

# below will output to both file(myoutput & mygraphs.pdf) and terminal
sink("output.log", append = TRUE, split = TRUE)
pdf("my_graphs.pdf")
source("../hello_world.r")

# now set back to only output to terminal
sink()
dev.off()
source("../hello_world.r")
