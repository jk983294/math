import numpy as np
import numpy.linalg as la
import math

weight = np.array([1.0, 1.0, 1.0, 1.0]).T
returns = np.array([0.1, 0.25, 0.35, 0.3]).T
alpha = returns
covariance = np.array([[1.2, 0.1, 0.0, 0.0], [0.1, 1.0, 0.0, 0.0],
                       [0.0, 0.0, 0.8, 0.0], [0.0, 0.0, 0.0, 0.8]])

print("unit exposure: ", np.dot(weight.T, returns))
std_dev = np.dot(np.dot(weight.T, covariance), weight)
print("std_dev: ", std_dev)

# 组合A
IR_A = math.sqrt(np.dot(np.dot(alpha.T, la.inv(covariance)), alpha))
total_risk = 1 / IR_A
residual_risk = 1 / IR_A
print("IR_A: ", IR_A, "total_risk: ", total_risk, "residual_risk: ",
      residual_risk)
