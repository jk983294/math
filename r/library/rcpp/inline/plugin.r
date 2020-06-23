# sudo apt-get install libgsl-dev

gslrng <- "
int seed = Rcpp::as<int>(par) ;
gsl_rng_env_setup();
gsl_rng *r = gsl_rng_alloc (gsl_rng_default);
gsl_rng_set (r, (unsigned long) seed);
double v = gsl_rng_get (r);
gsl_rng_free(r);
return Rcpp::wrap(v);
"

plug <- Rcpp:::Rcpp.plugin.maker(include.before = "#include <gsl/gsl_rng.h>", libs = paste("-L/usr/local/lib/R/site-library/Rcpp/libs -lRcpp -Wl,-rpath,/usr/local/lib/R/site-library/Rcpp/lib ", 
    "-L/usr/lib -lgsl -lgslcblas -lm", sep = ""))

registerPlugin("gslDemo", plug)

fun <- cxxfunction(signature(par = "numeric"), gslrng, plugin = "gslDemo")

fun(0)
fun(42)
