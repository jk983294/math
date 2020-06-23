library(Rcpp)
library(inline)

# x(t) = Ax(t-1) + u(t) parameter and error terms used throughout
a <- matrix(c(0.5, 0.1, 0.1, 0.5), nrow = 2)
u <- matrix(rnorm(10000), ncol = 2)

rSim <- function(coeff, errors) {
    simdata <- matrix(0, nrow(errors), ncol(errors))
    for (row in 2:nrow(errors)) {
        simdata[row, ] = coeff %*% simdata[(row - 1), ] + errors[row, ]
    }
    return(simdata)
}

rData <- rSim(a, u)

code <- "
    arma::mat coeff = Rcpp::as<arma::mat>(a);
    arma::mat errors = Rcpp::as<arma::mat>(u);
    int m = errors.n_rows;
    int n = errors.n_cols;
    arma::mat simdata(m,n);
    simdata.row(0) = arma::zeros<arma::mat>(1,n);
    for (int row=1; row<m; row++) {
        simdata.row(row) = simdata.row(row-1)*trans(coeff)
        + errors.row(row);
    }
    return Rcpp::wrap(simdata);
"

rcppSim <- cxxfunction(signature(a = "numeric", u = "numeric"), code, plugin = "RcppArmadillo")

rcppData <- rcppSim(a, u)

stopifnot(all.equal(rData, rcppData))


library(rbenchmark)
benchmark(rcppSim(a, u), rSim(a, u), columns = c("test", "replications", "elapsed", 
    "relative", "user.self", "sys.self"), order = "relative")
