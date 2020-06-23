require(Rcpp)
require(inline)

inc <- "
using namespace Rcpp;
double norm( double x, double y ) {
    return sqrt( x*x + y*y );
}
std::string hello() {
    return \"hello\";
}
int bar( int x) {
    return x*2;
}
double foo( int x, double y) {
    return x * y;
}
void bla( ) {
    Rprintf( \"hello \\n\" );
}
void bla1( int x) {
    Rprintf( \"hello (x = %d) \\n\", x );
}
void bla2( int x, double y) {
    Rprintf( \"hello (x = %d, y = %5.2f) \\n\", x, y );
}
RCPP_MODULE(mod) {
    function( \"norm\", &norm, List::create( _[\"x\"] = 0.0, _[\"y\"] = 0.0 ), \"Provides a simple vector norm\" );
    function(\"hello\", &hello);
    function(\"bar\", &bar );
    function(\"foo\", &foo, List::create( _[\"x\"], _[\"...\"] ), \"... means optional arguments\" );
    function(\"bla\", &bla );
    function(\"bla1\", &bla1 );
    function(\"bla2\", &bla2, List::create( _[\"x\"], _[\"y\"] = 0.0 ), \"x do not have default value\" );
}
"

fx <- cxxfunction(signature(), plugin = "Rcpp", include = inc)
mod <- Module("mod", getDynLib(fx))
show(mod$norm)
args(mod$norm)
mod$norm()  # use default value
mod$norm(y = 2)
mod$norm(3, 4)
mod$bar(2L)
mod$foo(2L, 10)
mod$hello()
mod$bla()
mod$bla1(2L)
mod$bla2(2L, 5)
