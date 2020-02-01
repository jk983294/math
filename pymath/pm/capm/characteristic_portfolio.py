import numpy as np
import numpy.linalg as la
import math

risk_free = 0.02
# weights = np.array([0.1, 0.2, 0.3, 0.3]).T
# betas = np.array([0.7, 0.9, 1.2, 1.2]).T
expected_returns = np.array([0.1, 0.2, 0.3, 0.3]).T
covariance = np.array([[1.2, 0.1, 0.0, 0.0], [0.1, 1.0, 0.0, 0.0],
                       [0.0, 0.0, 0.8, 0.0], [0.0, 0.0, 0.0, 0.8]])
capitals = np.array([0.5, 0.3, 0.1, 0.1]).T  # characteristics / attributes

# given characteristics, we calculate minimum risk
cp_weights = np.dot(la.inv(covariance), capitals) / np.dot(
    capitals.T, np.dot(la.inv(covariance), capitals))
cp_variance = 1 / np.dot(capitals.T, np.dot(la.inv(covariance), capitals))
cp_betas = np.dot(
    covariance, cp_weights) / cp_variance  # exactly = capitals characteristics
print('cp_weights', cp_weights)
print('total_weights', cp_weights.sum())
print('cp_variance', cp_variance)
print('cp_stddev', math.sqrt(cp_variance))
print('cp_betas', cp_betas)
