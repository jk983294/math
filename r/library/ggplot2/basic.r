library(ggplot2)

ggplot(
    data = BOD, # data
    mapping = aes(x = Time, y = demand) # mapping
) +
    geom_point(size = 5) + # geometric representation
    geom_line(color = "red") # another layer

# short format of above
ggplot(BOD, aes(Time, demand)) +
    geom_point(size = 5) +
    geom_line(color = "red")

ggplot(CO2, aes(x = conc, y = uptake, color = Treatment)) +
    geom_point(size = 4) +
    geom_smooth(method = lm, se = FALSE) +
    facet_wrap(~Type) + # group by Type
    labs(title = "CO2") +
    theme_bw()

ggplot(CO2, aes(Treatment, uptake)) +
    geom_boxplot() +
    geom_point(
        alpha = 0.5,
        aes(size = conc, color = Plant)
    ) # aes only apply to this point layer
