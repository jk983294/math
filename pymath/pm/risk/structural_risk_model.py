import numpy as np
import numpy.linalg as la
import math

# three factor, A industry, B industry and capital, 4 assets
# first column of factor_loading means asset 1 & 2 belong to industry A, 3 & 4 do not belong to industry A
factor_loading = np.array([[1.0, 0.0, 0.5], [1.0, 0.0, -0.5], [0.0, 1.0, 0.5],
                           [0.0, 1.0, -0.5]])
factor_return = np.array([0.1, 0.2, 0.05]).T
specific_return = np.array([0.05, 0.03, 0.02, 0.01]).T
factor_return_cov = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0],
                              [0.0, 0.0, 1.0]])
# specific_return_cov is diagonal matrix
specific_return_cov = np.array([[0.3, 0.0, 0.0, 0.0], [0.0, 0.4, 0.0, 0.0],
                                [0.0, 0.0, 0.5, 0.0], [0.0, 0.0, 0.0, 0.6]])

# return decomposition
expected_return = np.dot(factor_loading, factor_return) + specific_return
print('expected_return', expected_return)

# estimate asset return covariance matrix
estimated_return_cov = np.dot(np.dot(factor_loading, factor_return_cov),
                              factor_loading.T) + specific_return_cov
print('estimated_return_cov\n', estimated_return_cov)

# with return observation, estimate factor return
return_observation = np.copy(
    expected_return
)  # this is reality, we don't have specific_return subtracted
# return_observation = expected_return - specific_return # # this gives exact factor_return
src_inv = la.inv(specific_return_cov)
tmp = la.inv(np.dot(np.dot(factor_loading.T, src_inv), factor_loading))
estimated_factor_return = np.dot(
    np.dot(np.dot(tmp, factor_loading.T), src_inv), return_observation)
print('estimated_factor_return\n', estimated_factor_return)
"""
active risk
"""
bench_weight = np.array([0.2, 0.3, 0.2, 0.3]).T
portfolio_weight = np.array([0.1, 0.4, 0.3, 0.2]).T
bench_exposure = np.dot(factor_loading.T, bench_weight)
portfolio_exposure = np.dot(factor_loading.T, portfolio_weight)
delta_weight = portfolio_weight - bench_weight
active_portfolio_exposure = np.dot(factor_loading.T, delta_weight)
bench_variance = np.dot(np.dot(bench_weight.T, estimated_return_cov),
                        bench_weight)
portfolio_variance = np.dot(np.dot(portfolio_weight.T, estimated_return_cov),
                            portfolio_weight)
active_variance = np.dot(np.dot(delta_weight.T, estimated_return_cov),
                         delta_weight)
print('active_variance', active_variance)

b = np.dot(factor_return_cov, bench_exposure) / bench_variance
d = np.dot(specific_return_cov, bench_weight) / bench_variance
portfolio_beta = np.dot(portfolio_exposure.T, b) + np.dot(
    portfolio_weight.T, d)
print('portfolio_beta', portfolio_beta)
"""
attribution of risk
"""
# asset marginal contribution total risk
asset_mctr = np.dot(estimated_return_cov,
                    portfolio_weight) / math.sqrt(portfolio_variance)
print('asset marginal_contribution_total_variance', asset_mctr)
asset_mcar = np.dot(estimated_return_cov,
                    delta_weight) / math.sqrt(active_variance)
print('asset marginal_contribution_active_risk', asset_mcar)
asset_rmcar = np.multiply(asset_mcar, delta_weight)
print('asset relative marginal_contribution_active_risk', asset_rmcar)

# factor marginal contribution
factor_mcar = np.dot(factor_return_cov,
                     active_portfolio_exposure) / math.sqrt(active_variance)
print('factor marginal_contribution_active_risk', factor_mcar)
