set.seed(101)

biased_outcomes <- sample(c(0, 1), 1000, replace = TRUE, prob = c(0.4, 0.6))
prob_estimate <- sum(biased_outcomes)/length(biased_outcomes)  # Frequentist
