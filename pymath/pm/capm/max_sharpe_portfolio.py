import numpy as np
import numpy.linalg as la
import math

risk_free = 0.02
expected_returns = np.array([0.1, 0.2, 0.3, 0.3]).T
covariance = np.array([[1.2, 0.1, 0.0, 0.0], [0.1, 1.0, 0.0, 0.0],
                       [0.0, 0.0, 0.8, 0.0], [0.0, 0.0, 0.0, 0.8]])

# now expected_returns is characteristics, we calculate minimum risk max return
cp_variance = 1 / np.dot(expected_returns.T,
                         np.dot(la.inv(covariance), expected_returns))
cp_weights = np.dot(la.inv(covariance), expected_returns) * cp_variance
cp_betas = np.dot(
    covariance,
    cp_weights) / cp_variance  # exactly = expected_returns characteristics
sharpe = math.sqrt(1 / cp_variance)
print('cp_weights', cp_weights)
print('total_weights', cp_weights.sum())
print('cp_variance', cp_variance)
print('cp_stddev', math.sqrt(cp_variance))
print('cp_betas', cp_betas)
print('sharpe', sharpe)
