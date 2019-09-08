set.seed(1234)  # make results reproducible

## Exploratory factor analysis of ability.cov data
library(psych)
covariances <- ability.cov$cov
# convert covariances to correlations
correlations <- cov2cor(covariances)
correlations

# determine number of factors to extract
fa.parallel(correlations, n.obs = 112, fa = "both", n.iter = 100, main = "Scree plots with parallel analysis") # 2 factor


# Principal axis factoring without rotation
fa <- fa(correlations, nfactors = 2, rotate = "none", fm = "pa")
fa


# Factor extraction with orthogonal rotation
fa.varimax <- fa(correlations, nfactors = 2, rotate = "varimax", fm = "pa") # 正交旋转将人为地强制两个因子不相关
fa.varimax


# Factor extraction with oblique rotation
fa.promax <- fa(correlations, nfactors = 2, rotate = "promax", fm = "pa")
fa.promax

# calculate factor loading matrix
fsm <- function(oblique) {
    if (class(oblique)[2] == "fa" & is.null(oblique$Phi)) {
        warning("Object doesn't look like oblique EFA")
    } else {
        P <- unclass(oblique$loading)
        F <- P %*% oblique$Phi
        colnames(F) <- c("PA1", "PA2")
        return(F)
    }
}
fsm(fa.promax)

# plot factor solution
factor.plot(fa.promax, labels = rownames(fa.promax$loadings))
fa.diagram(fa.promax, simple = FALSE)

# factor scores
fa.promax$weights
