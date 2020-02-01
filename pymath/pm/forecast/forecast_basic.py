import numpy as np
import numpy.linalg as la
import math

# one asset k forecast
asset_return = 0.04
asset_std = 0.1
forecasts = np.array([0.05, 0.02, 0.03]).T
ICs = np.array([0.05, 0.06, 0.10]).T
forecast_corr = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
forecast_std = np.array([0.15, 0.06, 0.13])

forecast_std_mat = np.diag(forecast_std)
forecast_var = np.dot(np.dot(forecast_std_mat, forecast_corr),
                      forecast_std_mat)
forecast_cov = asset_std * np.dot(ICs.T, forecast_std_mat)
print("forecast_var: \n", forecast_var)
print("forecast_cov: \n", forecast_cov)

adjust_forecasts = forecasts - np.mean(forecasts)
differ_adjust = np.dot(np.dot(forecast_cov, la.inv(forecast_corr)),
                       adjust_forecasts)
print("differ_adjust: ", differ_adjust)

IC_portfolio = math.sqrt(np.dot(np.dot(ICs, la.inv(forecast_corr)), ICs.T))
print("IC_portfolio: ", IC_portfolio)
