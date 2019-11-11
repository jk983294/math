my_formula <- as.formula("y ~ x")
my_formula_without_intercept <- as.formula("y ~ x - 1")
class(my_formula)  # 'formula'
?formula
