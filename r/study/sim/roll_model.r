require(zoo)

trial <- 1000  # Number of trial
cost <- c()  # cost each trial
sd <- c()  # sd each trial

true.cost = 2
true.sd = 1

time = 1:trial

for (i in 1:trial) {
    # simulated Price Series
    epsilon = rnorm(time, sd = true.sd)
    prices = cumsum(epsilon)
    m_t = zoo(prices)
    a_t = m_t + true.cost
    b_t = m_t - true.cost
    
    # simulated trade prices
    q_t = sign(rnorm(time))
    p_t = m_t + (true.cost * q_t)
    
    # 1st difference of prices
    delta_p <- p_t - lag(p_t)
    # omit n.a. entry
    delta_p <- na.omit(delta_p)
    
    gamma_0 <- var(delta_p)
    gamma_1 <- cov(delta_p[1:length(delta_p) - 1], delta_p[2:length(delta_p)])
    sigma_2 <- gamma_0 + 2 * gamma_1
    
    if (gamma_1 > 0) {
        print("Error: Positive Autocovariance!")
    } else {
        cost <- append(cost, sqrt(-1 * gamma_1))
        sd <- append(sd, sigma_2)
    }
}

# Stimulated Cost Plot
plot(cost)
est.cost <- mean(cost)

plot(sd)
est.sd <- mean(sd)

# Final Result
cat("True cost and sd are", true.cost, " and ", true.sd)
cat("Estimated cost and sd are", est.cost, " and ", est.sd)
