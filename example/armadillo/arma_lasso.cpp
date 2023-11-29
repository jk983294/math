#include <math_stats.h>
#include <armadillo>
#include <cfloat>
#include <iostream>

using namespace std;

namespace mlpack {

// beta is the estimator
// yHat is the prediction from the current estimator

/**
 * An implementation of LARS, a stage-wise homotopy-based algorithm for
 * l1-regularized linear regression (LASSO) and l1+l2 regularized linear
 * regression (Elastic Net).
 *
 * Let \f$ X \f$ be a matrix where each row is a point and each column is a
 * dimension and let \f$ y \f$ be a vector of responses.
 *
 * The Elastic Net problem is to solve
 *
 * \f[ \min_{\beta} 0.5 || X \beta - y ||_2^2 + \lambda_1 || \beta ||_1 +
 *     0.5 \lambda_2 || \beta ||_2^2 \f]
 *
 * where \f$ \beta \f$ is the vector of regression coefficients.
 *
 * If \f$ \lambda_1 > 0 \f$ and \f$ \lambda_2 = 0 \f$, the problem is the LASSO.
 * If \f$ \lambda_1 > 0 \f$ and \f$ \lambda_2 > 0 \f$, the problem is the elastic net.
 * If \f$ \lambda_1 = 0 \f$ and \f$ \lambda_2 > 0 \f$, the problem is ridge regression.
 * If \f$ \lambda_1 = 0 \f$ and \f$ \lambda_2 = 0 \f$, the problem is unregularized linear regression.
 *
 * Note: This algorithm is not recommended for use (in terms of efficiency)
 * when \f$ \lambda_1 \f$ = 0.
 */
class LARS_arma {
public:
    /**
     * Set the parameters to LARS.  Both lambda1 and lambda2 default to 0.
     *
     * @param useCholesky Whether or not to use Cholesky decomposition when
     *    solving linear system (as opposed to using the full Gram matrix).
     * @param lambda1 Regularization parameter for l1-norm penalty.
     * @param lambda2 Regularization parameter for l2-norm penalty.
     * @param tolerance Run until the maximum correlation of elements in (X^T y)
     *     is less than this.
     * @param intercept If true, fit an intercept in the model.
     * @param normalize If true, normalize all features to have unit variance for
     *     training.
     * @param fitIntercept If true, fit an intercept in the model.
     * @param normalizeData If true, normalize all features to have unit variance
     * for training.
     */
    LARS_arma(bool useCholesky = false, double lambda1 = 0.0, double lambda2 = 0.0,
              double tolerance = 1e-16, bool fitIntercept = true, bool normalizeData = true);

    /**
     * Set the parameters to LARS, and pass in a precalculated Gram matrix.  Both
     * lambda1 and lambda2 default to 0.
     *
     * Note that the precalculated Gram matrix must match the settings of
     * `fitIntercept` and `normalizeData` (which both default to `true`): so, this
     * means that by default, the Gram matrix should be computed on mean-centered
     * data whose features are normalized to have unit variance.
     *
     * @param useCholesky Whether or not to use Cholesky decomposition when
     *    solving linear system (as opposed to using the full Gram matrix).
     * @param gramMatrix Gram matrix.
     * @param lambda1 Regularization parameter for l1-norm penalty.
     * @param lambda2 Regularization parameter for l2-norm penalty.
     * @param tolerance Run until the maximum correlation of elements in (X^T y)
     *     is less than this.
     * @param intercept If true, fit an intercept in the model.
     * @param normalize If true, normalize all features to have unit variance for
     *     training.
     * @param fitIntercept If true, fit an intercept in the model.
     * @param normalizeData If true, normalize all features to have unit variance
     * for training.
     */
    LARS_arma(const bool useCholesky, const arma::mat& gramMatrix, double lambda1 = 0.0,
              double lambda2 = 0.0, double tolerance = 1e-16, bool fitIntercept = true,
              bool normalizeData = true);

    /**
     * Set the parameters to LARS and run training. Both lambda1 and lambda2
     * are set by default to 0.
     *
     * @param data Input data.
     * @param responses A vector of targets.
     * @param transposeData Should be true if the input data is column-major and
     *     false otherwise.
     * @param useCholesky Whether or not to use Cholesky decomposition when
     *     solving linear system (as opposed to using the full Gram matrix).
     * @param lambda1 Regularization parameter for l1-norm penalty.
     * @param lambda2 Regularization parameter for l2-norm penalty.
     * @param tolerance Run until the maximum correlation of elements in (X^T y)
     *     is less than this.
     * @param fitIntercept If true, fit an intercept in the model.
     * @param normalizeData If true, normalize all features to have unit variance
     * for training.
     */
    LARS_arma(const arma::mat& data, const arma::rowvec& responses, bool transposeData = true,
              bool useCholesky = false, double lambda1 = 0.0, double lambda2 = 0.0,
              double tolerance = 1e-16, bool fitIntercept = true, bool normalizeData = true);

    /**
     * Set the parameters to LARS, pass in a precalculated Gram matrix, and run
     * training. Both lambda1 and lambda2 are set by default to 0.
     *
     * Note that the precalculated Gram matrix must match the settings of
     * `fitIntercept` and `normalizeData` (which both default to `true`): so, this
     * means that by default, the Gram matrix should be computed on mean-centered
     * data whose features are normalized to have unit variance.
     *
     * @param data Input data.
     * @param responses A vector of targets.
     * @param transposeData Should be true if the input data is column-major and
     *     false otherwise.
     * @param useCholesky Whether or not to use Cholesky decomposition when
     *     solving linear system (as opposed to using the full Gram matrix).
     * @param gramMatrix Gram matrix.
     * @param lambda1 Regularization parameter for l1-norm penalty.
     * @param lambda2 Regularization parameter for l2-norm penalty.
     * @param tolerance Run until the maximum correlation of elements in (X^T y)
     *     is less than this.
     * @param fitIntercept If true, fit an intercept in the model.
     * @param normalizeData If true, normalize all features to have unit variance
     * for training.
     */
    LARS_arma(const arma::mat& data, const arma::rowvec& responses, bool transposeData, bool useCholesky,
              const arma::mat& gramMatrix, double lambda1 = 0.0, double lambda2 = 0.0,
              double tolerance = 1e-16, bool fitIntercept = true, bool normalizeData = true);

    /**
     * Run LARS.  The input matrix (like all mlpack matrices) should be
     * column-major -- each column is an observation and each row is a dimension.
     * However, because LARS is more efficient on a row-major matrix, this method
     * will (internally) transpose the matrix.  If this transposition is not
     * necessary (i.e., you want to pass in a row-major matrix), pass 'false' for
     * the transposeData parameter.
     *
     * @param data Column-major input data (or row-major input data if rowMajor =
     *     true).
     * @param responses A vector of targets.
     * @param beta Vector to store the solution (the coefficients) in.
     * @param transposeData Set to false if the data is row-major.
     * @return minimum cost error(||y-beta*X||2 is used to calculate error).
     */
    double Train(const arma::mat& data, const arma::rowvec& responses, arma::vec& beta,
                 bool transposeData = true);

    /**
     * Run LARS.  The input matrix (like all mlpack matrices) should be
     * column-major -- each column is an observation and each row is a dimension.
     * However, because LARS is more efficient on a row-major matrix, this method
     * will (internally) transpose the matrix.  If this transposition is not
     * necessary (i.e., you want to pass in a row-major matrix), pass 'false' for
     * the transposeData parameter.
     *
     * @param data Input data.
     * @param responses A vector of targets.
     * @param transposeData Should be true if the input data is column-major and
     *     false otherwise.
     * @return minimum cost error(||y-beta*X||2 is used to calculate error).
     */
    double Train(const arma::mat& data, const arma::rowvec& responses, bool transposeData = true);

    /**
     * Predict y_i for each data point in the given data matrix using the
     * currently-trained LARS model.
     *
     * @param points The data points to regress on.
     * @param predictions y, which will contained calculated values on completion.
     * @param rowMajor Should be true if the data points matrix is row-major and
     *     false otherwise.
     */
    void Predict(const arma::mat& points, arma::rowvec& predictions, bool rowMajor = false) const;

    double Lambda1() const { return lambda1; }
    double& Lambda1() { return lambda1; }
    double Lambda2() const { return lambda2; }
    double& Lambda2() { return lambda2; }
    bool UseCholesky() const { return useCholesky; }
    bool& UseCholesky() { return useCholesky; }
    double Tolerance() const { return tolerance; }
    double& Tolerance() { return tolerance; }
    bool FitIntercept() const { return fitIntercept; }
    //! Modify whether or not to fit an intercept.
    void FitIntercept(const bool newFitIntercept) {
        // If we are storing a Gram matrix internally, but now will be normalizing
        // data, then the Gram matrix we have computed is incorrect and needs to be
        // recomputed.
        if (fitIntercept != newFitIntercept) {
            if (matGram != &matGramInternal) {
                throw std::invalid_argument(
                    "LARS::FitIntercept(): cannot change value "
                    "when an external Gram matrix was specified!");
            }

            fitIntercept = newFitIntercept;
            matGramInternal.clear();
        }
    }

    bool NormalizeData() const { return normalizeData; }
    //! Modify whether or not to normalize data during training.
    void NormalizeData(const bool newNormalizeData) {
        // If we are storing a Gram matrix internally, but now will be normalizing
        // data, then the Gram matrix we have computed is incorrect and needs to be
        // recomputed.
        if (normalizeData != newNormalizeData) {
            if (matGram != &matGramInternal) {
                throw std::invalid_argument(
                    "LARS::NormalizeData(): cannot change value "
                    "when an external Gram matrix was specified!");
            }

            normalizeData = newNormalizeData;
            matGramInternal.clear();
        }
    }

    //! Access the set of active dimensions.
    const std::vector<size_t>& ActiveSet() const { return activeSet; }

    //! Access the set of coefficients after each iteration; the solution is the
    //! last element.
    const std::vector<arma::vec>& BetaPath() const { return betaPath; }

    //! Access the solution coefficients
    const arma::vec& Beta() const { return betaPath.back(); }

    //! Access the set of values for lambda1 after each iteration; the solution is
    //! the last element.
    const std::vector<double>& LambdaPath() const { return lambdaPath; }

    //! Return the intercept (if fitted, otherwise 0).
    double Intercept() const { return interceptPath.back(); }

    //! Return the intercept path (the intercept for every model).
    const std::vector<double>& InterceptPath() const { return interceptPath; }

    //! Access the upper triangular cholesky factor.
    const arma::mat& MatUtriCholFactor() const { return matUtriCholFactor; }

    /**
     * Compute cost error of the given data matrix using the
     * currently-trained LARS model. Only ||y-beta*X||2 is used to calculate
     * cost error.
     *
     * @param matX Column-major input data (or row-major input data if rowMajor =
     *     true).
     * @param y responses A vector of targets.
     * @param rowMajor Should be true if the data points matrix is row-major and
     *   false otherwise.
     * @return The minimum cost error.
     */
    double ComputeError(const arma::mat& matX, const arma::rowvec& y, const bool rowMajor = false);

private:
    //! Gram matrix.
    arma::mat matGramInternal;

    //! Pointer to the Gram matrix we will use.
    const arma::mat* matGram;

    //! Upper triangular cholesky factor; initially 0x0 matrix.
    arma::mat matUtriCholFactor;

    //! Whether or not to use Cholesky decomposition when solving linear system.
    bool useCholesky;

    //! True if this is the LASSO problem.
    bool lasso;
    //! Regularization parameter for l1 penalty.
    double lambda1;

    //! True if this is the elastic net problem.
    bool elasticNet;
    //! Regularization parameter for l2 penalty.
    double lambda2;

    //! Tolerance for main loop.
    double tolerance;

    //! Whether or not to fit an intercept.
    bool fitIntercept;

    //! Whether or not to normalize data during training (i.e. make each feature
    //! have unit variance for training).
    bool normalizeData;

    //! Solution path.
    std::vector<arma::vec> betaPath;

    //! Value of lambda_1 for each solution in solution path.
    std::vector<double> lambdaPath;

    //! Intercept (only if fitIntercept is true).
    std::vector<double> interceptPath;

    //! Active set of dimensions.
    std::vector<size_t> activeSet;

    //! Active set membership indicator (for each dimension).
    std::vector<bool> isActive;

    // Set of variables that are ignored (if any).

    //! Set of ignored variables (for dimensions in span{active set dimensions}).
    std::vector<size_t> ignoreSet;

    //! Membership indicator for set of ignored variables.
    std::vector<bool> isIgnored;

    /**
     * Remove activeVarInd'th element from active set.
     *
     * @param activeVarInd Index of element to remove from active set.
     */
    void Deactivate(const size_t activeVarInd);

    /**
     * Add dimension varInd to active set.
     *
     * @param varInd Dimension to add to active set.
     */
    void Activate(const size_t varInd);

    /**
     * Add dimension varInd to ignores set (never removed).
     *
     * @param varInd Dimension to add to ignores set.
     */
    void Ignore(const size_t varInd);

    // compute "equiangular" direction in output space
    void ComputeYHatDirection(const arma::mat& matX, const arma::vec& betaDirection, arma::vec& yHatDirection);

    // interpolate to compute last solution vector
    void InterpolateBeta();

    void CholeskyInsert(const arma::vec& newX, const arma::mat& X);

    void CholeskyInsert(double sqNormNewX, const arma::vec& newGramCol);

    void GivensRotate(const arma::vec::fixed<2>& x, arma::vec::fixed<2>& rotatedX, arma::mat& G);

    void CholeskyDelete(const size_t colToKill);
};

inline LARS_arma::LARS_arma(bool useCholesky, double lambda1, double lambda2, double tolerance,
                            bool fitIntercept, bool normalizeData)
    : matGram(&matGramInternal),
      useCholesky(useCholesky),
      lasso((lambda1 != 0)),
      lambda1(lambda1),
      elasticNet((lambda1 != 0) && (lambda2 != 0)),
      lambda2(lambda2),
      tolerance(tolerance),
      fitIntercept(fitIntercept),
      normalizeData(normalizeData) { /* Nothing left to do. */
}

inline LARS_arma::LARS_arma(bool useCholesky, const arma::mat& gramMatrix, double lambda1,
                            double lambda2, double tolerance, bool fitIntercept,
                            bool normalizeData)
    : matGram(&gramMatrix),
      useCholesky(useCholesky),
      lasso((lambda1 != 0)),
      lambda1(lambda1),
      elasticNet((lambda1 != 0) && (lambda2 != 0)),
      lambda2(lambda2),
      tolerance(tolerance),
      fitIntercept(fitIntercept),
      normalizeData(normalizeData) { /* Nothing left to do */
}

inline LARS_arma::LARS_arma(const arma::mat& data, const arma::rowvec& responses, bool transposeData,
                            bool useCholesky, double lambda1, double lambda2, double tolerance,
                            bool fitIntercept, bool normalizeData)
    : LARS_arma(useCholesky, lambda1, lambda2, tolerance, fitIntercept, normalizeData) {
    Train(data, responses, transposeData);
}

inline LARS_arma::LARS_arma(const arma::mat& data, const arma::rowvec& responses, bool transposeData,
                            bool useCholesky, const arma::mat& gramMatrix, double lambda1,
                            double lambda2, double tolerance, bool fitIntercept,
                            bool normalizeData)
    : LARS_arma(useCholesky, gramMatrix, lambda1, lambda2, tolerance, fitIntercept, normalizeData) {
    Train(data, responses, transposeData);
}

inline double LARS_arma::Train(const arma::mat& matX, const arma::rowvec& y, arma::vec& beta,
                               bool transposeData) {
    // Clear any previous solution information.
    betaPath.clear();
    lambdaPath.clear();
    activeSet.clear();
    isActive.clear();
    ignoreSet.clear();
    isIgnored.clear();
    matUtriCholFactor.reset();

    // Update values in case lambda1 or lambda2 changed.
    lasso = (lambda1 != 0);
    elasticNet = (lambda1 != 0 && lambda2 != 0);

    // This matrix may end up holding the transpose -- if necessary.
    arma::mat dataTrans;
    // This vector may hold zero-centered responses, if necessary.
    arma::rowvec yCentered;

    // dataRef is row-major.  We can reuse the given matX, but only if we don't
    // need to do any transformations to it.
    const arma::mat& dataRef = (transposeData || fitIntercept || normalizeData) ? dataTrans : matX;
    const arma::rowvec& yRef = (fitIntercept) ? yCentered : y;

    arma::vec offsetX;     // used only if fitting an intercept
    double offsetY = 0.0;  // used only if fitting an intercept
    arma::vec stdX;        // used only if normalizing

    if (transposeData) {
        if (fitIntercept) {
            offsetX = arma::mean(matX, 1);
            dataTrans = (matX.each_col() - offsetX).t();
        }

        if (normalizeData) {
            stdX = arma::stddev(matX, 0, 1);
            stdX.replace(0.0, 1.0);  // Make sure we don't divide by 0!

            // Check if we have already done the transposition.
            if (!fitIntercept)
                dataTrans = (matX.each_col() / stdX).t();
            else
                dataTrans.each_row() /= stdX.t();
        }

        // Make sure we convert the data to row-major format, even if no
        // transformations were needed.
        if (!fitIntercept && !normalizeData) dataTrans = matX.t();
    } else {
        // We don't need to transpose the data---it's already in row-major form.
        if (fitIntercept) {
            offsetX = arma::mean(matX, 0).t();
            dataTrans = (matX.each_row() - offsetX.t());
        }

        if (normalizeData) {
            stdX = arma::stddev(matX, 0, 0).t();
            stdX.replace(0.0, 1.0);  // Make sure we don't divide by 0!

            // Check if we have already populated `dataTrans`.
            if (!fitIntercept)
                dataTrans = (matX.each_row() / stdX.t());
            else
                dataTrans.each_row() /= stdX.t();
        }

        // If we are not fitting an intercept and we are not normalizing the data,
        // dataTrans already points to matX so we don't need to do anything.
    }

    if (fitIntercept) {
        offsetY = arma::mean(y);
        yCentered = y - offsetY;
    }

    // Compute X' * y.
    arma::vec vecXTy = trans(yRef * dataRef);

    // Set up active set variables.  In the beginning, the active set has size 0
    // (all dimensions are inactive).
    isActive.resize(dataRef.n_cols, false);

    // Set up ignores set variables. Initialized empty.
    isIgnored.resize(dataRef.n_cols, false);

    // Initialize yHat and beta.
    beta = arma::zeros(dataRef.n_cols);
    arma::vec yHat = arma::zeros(dataRef.n_rows);
    arma::vec yHatDirection(dataRef.n_rows);

    bool lassocond = false;

    // Compute the initial maximum correlation among all dimensions.
    arma::vec corr = vecXTy;
    double maxCorr = 0;
    size_t changeInd = 0;
    size_t lassocondInd = dataRef.n_cols;
    for (size_t i = 0; i < vecXTy.n_elem; ++i) {
        if (fabs(corr(i)) > maxCorr) {
            maxCorr = fabs(corr(i));
            changeInd = i;
        }
    }

    betaPath.push_back(beta);
    lambdaPath.push_back(maxCorr);

    // If the maximum correlation is too small, there is no reason to continue.
    if (maxCorr < lambda1) {
        lambdaPath[0] = lambda1;

        if (fitIntercept)
            interceptPath.push_back(offsetY - arma::dot(offsetX, betaPath[0]));
        else
            interceptPath.push_back(0.0);

        return maxCorr;
    }

    // Compute the Gram matrix.  If this is the elastic net problem, we will add
    // lambda2 * I_n to the matrix.
    if (matGram->n_elem != dataRef.n_cols * dataRef.n_cols) {
        // In this case, matGram should reference matGramInternal.
        matGramInternal = trans(dataRef) * dataRef;

        if (elasticNet && !useCholesky) matGramInternal += lambda2 * arma::eye(dataRef.n_cols, dataRef.n_cols);
    }

    // Main loop.
    while (((activeSet.size() + ignoreSet.size()) < dataRef.n_cols) && (maxCorr > tolerance)) {
        // Compute the maximum correlation among inactive dimensions.
        maxCorr = 0;
        double maxActiveCorr = 0;
        double minActiveCorr = DBL_MAX;
        for (size_t i = 0; i < dataRef.n_cols; ++i) {
            if ((!isActive[i]) && (!isIgnored[i]) && (fabs(corr(i)) > maxCorr)) {
                maxCorr = fabs(corr(i));
                changeInd = i;
            } else if (isActive[i] && (matGram != &matGramInternal)) {
                // Here we will do a sanity check: if the correlation of any dimension
                // is not the maximum correlation, then the user has probably passed a
                // Gram matrix whose properties do not match the value of fitIntercept
                // and normalizeData.
                if (fabs(corr(i)) > maxActiveCorr) maxActiveCorr = fabs(corr(i));
                if (fabs(corr(i)) < minActiveCorr) minActiveCorr = fabs(corr(i));
            }
        }

        // If the maximum correlation is sufficiently small, don't add this
        // variable; terminate early.
        if (maxCorr < tolerance) break;

        if ((matGram != &matGramInternal) &&
            ((maxActiveCorr - minActiveCorr) > 100 * std::numeric_limits<double>::epsilon())) {
            // Construct the error message to match the user's settings.
            std::ostringstream oss;
            oss << "LARS::Train(): correlation conditions violated; check that your "
                << "given Gram matrix is properly computed on ";
            if (fitIntercept)
                oss << "mean-centered ";
            else
                oss << "non-mean-centered ";
            if (normalizeData)
                oss << "unit-variance (normalized) ";
            else
                oss << "non-normalized ";
            oss << "data!";
            oss << std::endl;
            throw std::runtime_error(oss.str());
        }

        // Add the variable to the active set and update the Gram matrix as
        // necessary.
        if (!lassocond) {
            if (useCholesky) {
                // vec newGramCol = vec(activeSet.size());
                // for (size_t i = 0; i < activeSet.size(); ++i)
                // {
                //   newGramCol[i] = dot(matX.col(activeSet[i]), matX.col(changeInd));
                // }
                // This is equivalent to the above 5 lines.
                arma::vec newGramCol =
                    matGram->elem(changeInd * dataRef.n_cols + arma::conv_to<arma::uvec>::from(activeSet));

                CholeskyInsert((*matGram)(changeInd, changeInd), newGramCol);
            }
            Activate(changeInd);
        }

        // Compute signs of correlations.
        arma::vec s = arma::vec(activeSet.size());
        for (size_t i = 0; i < activeSet.size(); ++i) s(i) = corr(activeSet[i]) / fabs(corr(activeSet[i]));

        // Compute the "equiangular" direction in parameter space (betaDirection).
        // We use quotes because in the case of non-unit norm variables, this need
        // not be equiangular.
        arma::vec unnormalizedBetaDirection;
        double normalization;
        arma::vec betaDirection;
        if (useCholesky) {
            // Check for singularity.
            const double lastUtriElement =
                matUtriCholFactor(matUtriCholFactor.n_cols - 1, matUtriCholFactor.n_rows - 1);
            if (std::abs(lastUtriElement) > tolerance) {
                // Ok, no singularity.
                /**
                 * Note that:
                 * R^T R % S^T % S = (R % S)^T (R % S)
                 * Now, for 1 the ones vector:
                 * inv( (R % S)^T (R % S) ) 1
                 *    = inv(R % S) inv((R % S)^T) 1
                 *    = inv(R % S) Solve((R % S)^T, 1)
                 *    = inv(R % S) Solve(R^T, s)
                 *    = Solve(R % S, Solve(R^T, s)
                 *    = s % Solve(R, Solve(R^T, s))
                 */
                unnormalizedBetaDirection =
                    solve(trimatu(matUtriCholFactor), solve(trimatl(trans(matUtriCholFactor)), s));

                normalization = 1.0 / sqrt(dot(s, unnormalizedBetaDirection));
                betaDirection = normalization * unnormalizedBetaDirection;
            } else {
                // Singularity, so remove variable from active set, add to ignores set,
                // and look for new variable to add.
                printf("Encountered singularity when adding variable %zu to active set; permanently removing.\n",
                       changeInd);
                Deactivate(activeSet.size() - 1);
                Ignore(changeInd);
                CholeskyDelete(matUtriCholFactor.n_rows - 1);

                // Note that although we are now ignoring this variable, we may still
                // need to take a step with the previous beta direction towards the next
                // variable we will add.
                s = s.subvec(0, activeSet.size() - 1);  // Drop last element.
                unnormalizedBetaDirection =
                    solve(trimatu(matUtriCholFactor), solve(trimatl(trans(matUtriCholFactor)), s));

                normalization = 1.0 / sqrt(dot(s, unnormalizedBetaDirection));
                betaDirection = normalization * unnormalizedBetaDirection;
            }
        } else {
            arma::mat matGramActive = arma::mat(activeSet.size(), activeSet.size());
            for (size_t i = 0; i < activeSet.size(); ++i)
                for (size_t j = 0; j < activeSet.size(); ++j)
                    matGramActive(i, j) = (*matGram)(activeSet[i], activeSet[j]);

            // Check for singularity.
            arma::mat matS = s * arma::ones<arma::mat>(1, activeSet.size());
            const bool solvedOk = solve(unnormalizedBetaDirection, matGramActive % trans(matS) % matS,
                                        arma::ones<arma::mat>(activeSet.size(), 1));
            if (solvedOk) {
                // Ok, no singularity.
                normalization = 1.0 / sqrt(sum(unnormalizedBetaDirection));
                betaDirection = normalization * unnormalizedBetaDirection % s;
            } else {
                // Singularity, so remove variable from active set, add to ignores set,
                // and look for new variable to add.
                Deactivate(activeSet.size() - 1);
                Ignore(changeInd);
                printf("Encountered singularity when adding variable %zu to active set; permanently removing.\n",
                       changeInd);

                // Note that although we are now ignoring this variable, we may still
                // need to take a step with the previous beta direction towards the next
                // variable we will add.
                s = s.subvec(0, activeSet.size() - 1);  // Drop last element.
                matS = s * arma::ones<arma::mat>(1, activeSet.size());
                // This worked last iteration, so there can't be a singularity.
                solve(unnormalizedBetaDirection, matGramActive % trans(matS) % matS,
                      arma::ones<arma::mat>(activeSet.size(), 1));
                normalization = 1.0 / sqrt(sum(unnormalizedBetaDirection));
                betaDirection = normalization * unnormalizedBetaDirection % s;
            }
        }

        // compute "equiangular" direction in output space
        ComputeYHatDirection(dataRef, betaDirection, yHatDirection);

        double gamma = maxCorr / normalization;

        // If not all variables are active.
        if ((activeSet.size() + ignoreSet.size()) < dataRef.n_cols) {
            // Compute correlations with direction.
            for (size_t ind = 0; ind < dataRef.n_cols; ind++) {
                if (isActive[ind] || isIgnored[ind]) continue;

                const double dirCorr = dot(dataRef.col(ind), yHatDirection);
                const double val1 = (maxCorr - corr(ind)) / (normalization - dirCorr);
                const double val2 = (maxCorr + corr(ind)) / (normalization + dirCorr);

                // If we kicked out a feature due to the LASSO modification last
                // iteration, then we do not allow relaxation of the step size to 0 for
                // that feature in this iteration.
                if (lassocond && (ind == lassocondInd)) {
                    if ((val1 > 0.0) && (val1 < gamma)) gamma = val1;
                    if ((val2 > 0.0) && (val2 < gamma)) gamma = val2;
                } else {
                    if ((val1 >= 0.0) && (val1 < gamma)) gamma = val1;
                    if ((val2 >= 0.0) && (val2 < gamma)) gamma = val2;
                }
            }
        }

        // Bound gamma according to LASSO.
        if (lasso) {
            lassocond = false;
            lassocondInd = dataRef.n_cols;
            double lassoboundOnGamma = DBL_MAX;
            size_t activeIndToKickOut = -1;

            for (size_t i = 0; i < activeSet.size(); ++i) {
                double val = -beta(activeSet[i]) / betaDirection(i);
                if ((val > 0) && (val < lassoboundOnGamma)) {
                    lassoboundOnGamma = val;
                    activeIndToKickOut = i;
                }
            }

            if (lassoboundOnGamma < gamma) {
                gamma = lassoboundOnGamma;
                lassocond = true;
                changeInd = activeIndToKickOut;
                lassocondInd = activeSet[changeInd];
            }
        }

        // Update the prediction.
        yHat += gamma * yHatDirection;

        // Update the estimator.
        for (size_t i = 0; i < activeSet.size(); ++i) {
            beta(activeSet[i]) += gamma * betaDirection(i);
        }

        // Sanity check to make sure the kicked out dimension is actually zero.
        if (lassocond) {
            if (beta(activeSet[changeInd]) != 0) {
                beta(activeSet[changeInd]) = 0;
            }
        }

        betaPath.push_back(beta);

        if (lassocond) {
            // Index is in position changeInd in activeSet.
            if (useCholesky) CholeskyDelete(changeInd);

            Deactivate(changeInd);
        }

        corr = vecXTy - trans(dataRef) * yHat;
        if (elasticNet) corr -= lambda2 * beta;

        double curLambda = 0;
        for (size_t i = 0; i < activeSet.size(); ++i) curLambda += fabs(corr(activeSet[i]));

        curLambda /= ((double)activeSet.size());

        lambdaPath.push_back(curLambda);

        // Time to stop for LASSO?
        if (lasso) {
            if (curLambda <= lambda1) {
                InterpolateBeta();
                break;
            }
        }
    }

    // Perform un-scaling of learned beta, if needed, to account for
    // normalization.
    if (normalizeData) {
        for (size_t i = 0; i < betaPath.size(); ++i) {
            betaPath[i] /= stdX;
        }
    }

    // Set the intercept values.  This is needed (for paranoia reasons) even if an
    // intercept isn't fit, in case a user changes `fitIntercept` after training.
    // If an intercept wasn't fit, we set them all to zero.
    if (fitIntercept) {
        interceptPath.clear();
        for (size_t i = 0; i < betaPath.size(); ++i) interceptPath.push_back(offsetY - arma::dot(offsetX, betaPath[i]));
    } else {
        interceptPath.clear();
        interceptPath.resize(betaPath.size(), 0.0);
    }

    // Unfortunate copy...
    beta = betaPath.back();

    return ComputeError(matX, y, !transposeData);
}

inline double LARS_arma::Train(const arma::mat& data, const arma::rowvec& responses, bool transposeData) {
    arma::vec beta;
    return Train(data, responses, beta, transposeData);
}

inline void LARS_arma::Predict(const arma::mat& points, arma::rowvec& predictions, bool rowMajor) const {
    // We really only need to store beta internally...
    if (rowMajor && !fitIntercept)
        predictions = trans(points * betaPath.back());
    else if (rowMajor)
        predictions = trans(points * betaPath.back()) + interceptPath.back();
    else if (fitIntercept)
        predictions = betaPath.back().t() * points + interceptPath.back();
    else
        predictions = betaPath.back().t() * points;
}

// Private functions.
inline void LARS_arma::Deactivate(const size_t activeVarInd) {
    isActive[activeSet[activeVarInd]] = false;
    activeSet.erase(activeSet.begin() + activeVarInd);
}

inline void LARS_arma::Activate(const size_t varInd) {
    isActive[varInd] = true;
    activeSet.push_back(varInd);
}

inline void LARS_arma::Ignore(const size_t varInd) {
    isIgnored[varInd] = true;
    ignoreSet.push_back(varInd);
}

inline void LARS_arma::ComputeYHatDirection(const arma::mat& matX, const arma::vec& betaDirection,
                                            arma::vec& yHatDirection) {
    yHatDirection.fill(0);
    for (size_t i = 0; i < activeSet.size(); ++i) yHatDirection += betaDirection(i) * matX.col(activeSet[i]);
}

inline void LARS_arma::InterpolateBeta() {
    const size_t pathLength = betaPath.size();

    // interpolate beta and stop
    double ultimateLambda = lambdaPath[pathLength - 1];
    double penultimateLambda = lambdaPath[pathLength - 2];
    double interp = (penultimateLambda - lambda1) / (penultimateLambda - ultimateLambda);

    betaPath[pathLength - 1] = (1 - interp) * (betaPath[pathLength - 2]) + interp * betaPath[pathLength - 1];

    lambdaPath[pathLength - 1] = lambda1;
}

inline void LARS_arma::CholeskyInsert(const arma::vec& newX, const arma::mat& X) {
    if (matUtriCholFactor.n_rows == 0) {
        matUtriCholFactor = arma::mat(1, 1);

        if (elasticNet)
            matUtriCholFactor(0, 0) = sqrt(dot(newX, newX) + lambda2);
        else
            matUtriCholFactor(0, 0) = norm(newX, 2);
    } else {
        arma::vec newGramCol = trans(X) * newX;
        CholeskyInsert(dot(newX, newX), newGramCol);
    }
}

inline void LARS_arma::CholeskyInsert(double sqNormNewX, const arma::vec& newGramCol) {
    int n = matUtriCholFactor.n_rows;

    if (n == 0) {
        matUtriCholFactor = arma::mat(1, 1);

        if (elasticNet)
            matUtriCholFactor(0, 0) = sqrt(sqNormNewX + lambda2);
        else
            matUtriCholFactor(0, 0) = sqrt(sqNormNewX);
    } else {
        arma::mat matNewR = arma::mat(n + 1, n + 1);

        if (elasticNet) sqNormNewX += lambda2;

        arma::vec matUtriCholFactork = solve(trimatl(trans(matUtriCholFactor)), newGramCol);

        matNewR(arma::span(0, n - 1), arma::span(0, n - 1)) = matUtriCholFactor;
        matNewR(arma::span(0, n - 1), n) = matUtriCholFactork;
        matNewR(n, arma::span(0, n - 1)).fill(0.0);
        matNewR(n, n) = sqrt(sqNormNewX - dot(matUtriCholFactork, matUtriCholFactork));

        matUtriCholFactor = matNewR;
    }
}

inline void LARS_arma::GivensRotate(const arma::vec::fixed<2>& x, arma::vec::fixed<2>& rotatedX, arma::mat& matG) {
    if (x(1) == 0) {
        matG = arma::eye(2, 2);
        rotatedX = x;
    } else {
        double r = norm(x, 2);
        matG = arma::mat(2, 2);

        double scaledX1 = x(0) / r;
        double scaledX2 = x(1) / r;

        matG(0, 0) = scaledX1;
        matG(1, 0) = -scaledX2;
        matG(0, 1) = scaledX2;
        matG(1, 1) = scaledX1;

        rotatedX = arma::vec(2);
        rotatedX(0) = r;
        rotatedX(1) = 0;
    }
}

inline void LARS_arma::CholeskyDelete(const size_t colToKill) {
    size_t n = matUtriCholFactor.n_rows;

    if (colToKill == (n - 1)) {
        matUtriCholFactor = matUtriCholFactor(arma::span(0, n - 2), arma::span(0, n - 2));
    } else {
        matUtriCholFactor.shed_col(colToKill);  // remove column colToKill
        n--;

        for (size_t k = colToKill; k < n; ++k) {
            arma::mat matG;
            arma::vec::fixed<2> rotatedVec;
            GivensRotate(matUtriCholFactor(arma::span(k, k + 1), k), rotatedVec, matG);
            matUtriCholFactor(arma::span(k, k + 1), k) = rotatedVec;
            if (k < n - 1) {
                matUtriCholFactor(arma::span(k, k + 1), arma::span(k + 1, n - 1)) =
                    matG * matUtriCholFactor(arma::span(k, k + 1), arma::span(k + 1, n - 1));
            }
        }

        matUtriCholFactor.shed_row(n);
    }
}

inline double LARS_arma::ComputeError(const arma::mat& matX, const arma::rowvec& y, const bool rowMajor) {
    if (rowMajor) {
        return arma::accu(arma::pow(y - trans(matX * betaPath.back()), 2.0));
    }

    else {
        return arma::accu(arma::pow(y - betaPath.back().t() * matX, 2.0));
    }
}
}  // namespace mlpack

int main(int argc, char* argv[]) {
    double lambda1 = 1;
    double lambda2 = 0;
    bool useCholesky = true;
    bool noIntercept = true;
    bool noNormalize = true;

    auto* lars = new mlpack::LARS_arma(useCholesky, lambda1, lambda2);
    lars->FitIntercept(!noIntercept);
    lars->NormalizeData(!noNormalize);

    std::random_device rd;
    std::mt19937 generator(rd());
    std::normal_distribution<double> nd(0., 0.001);

    // y = 2 * x0 + 3 * x1  - 1 + delta
    size_t n = 1000;
    arma::mat matX(n * n, 3);
    arma::rowvec y(n * n);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            double x0 = i, x1 = j;
            double y_val = 2 * x0 + 3 * x1 + nd(generator);
            y(i * n + j) = y_val;
            matX(i * n + j, 0) = x0;
            matX(i * n + j, 1) = x1;
            matX(i * n + j, 2) = nd(generator) * 100000;
        }
    }

    if (y.n_rows > 1) cout << "Only one column or row allowed in responses file!" << endl;
    if (y.n_elem != matX.n_rows) cout << "Number of responses must be equal to number of rows of X!" << endl;

    arma::vec beta;
    lars->Train(matX, y, beta, false /* do not transpose */);

    arma::rowvec predictions;
    lars->Predict(matX.t(), predictions, false);
    cout << beta << endl;

    double pcor = ornate::corr(predictions.memptr(), y.memptr(), n * n);
    printf("corr %f\n", pcor);
}
