require( Rcpp )
require( inline )

inc <- '
using namespace Rcpp;
class Uniform {
public :
    Uniform(double min_, double max_) : min(min_), max(max_) {}
    NumericVector draw(int n) const {
        RNGScope scope;
        return runif( n, min, max );
    }
    double min, max;
};

double uniformRange( Uniform* w) {
    return w->max - w->min;
}

class Foo {
public :
    Foo(double x_, double y_, double z_ ): x(x_), y(y_), z(z_) {}
    double x;
    double y;

    double get_z() { return z; }
    void set_z( double z_ ) { z = z_; }
private :
    double z;
};

RCPP_MODULE(unif_module) {
    class_<Uniform>( "Uniform" )
        .constructor<double,double>("sets the min and max value of the distribution")
        .field( "min", &Uniform::min )
        .field( "max", &Uniform::max )
        .method( "draw", &Uniform::draw )
        .method( "range", &uniformRange )
        ;

    class_<Foo>( "Foo" )
        .constructor<double,double,double>()
        .field( "x", &Foo::x, "documentation for x" )
        .field_readonly( "y", &Foo::y, "documentation for y" )
        .property( "z", &Foo::get_z, &Foo::set_z, "Documentation for z" )
        ;
}
'

fx_unif <- cxxfunction(signature(), plugin="Rcpp", include=inc)
unif_module <- Module( "unif_module", getDynLib(fx_unif ) )
Uniform <- unif_module$Uniform
u <- new( Uniform, 0, 10 )
u$draw( 10L )
u$range()
u$max <- 1
u$range()
u$draw( 10 )

inc_std_vec <- '

typedef std::vector<double> vec; // convenience typedef
void vec_assign( vec* obj, Rcpp::NumericVector data ) { // helpers
    obj->assign( data.begin(), data.end() );
}

void vec_insert( vec* obj, int position, Rcpp::NumericVector data) {
    vec::iterator it = obj->begin() + position;
    obj->insert( it, data.begin(), data.end() );
}

Rcpp::NumericVector vec_asR( vec* obj ) { return Rcpp::wrap( *obj ); }

void vec_set( vec* obj, int i, double value ) { obj->at( i ) = value; }

RCPP_MODULE(mod_vec) {
    using namespace Rcpp;

    // we expose the class std::vector < double > as "vec" on the R side
    class_<vec>( "vec")
        // exposing constructors
        .constructor()
        .constructor<int>()
        // exposing member functions
        .method( "size", &vec::size)
        .method( "max_size", &vec::max_size)
        .method( "resize", &vec::resize)
        .method( "capacity", &vec::capacity)
        .method( "empty", &vec::empty)
        .method( "reserve", &vec::reserve)
        .method( "push_back", &vec::push_back )
        .method( "pop_back", &vec::pop_back )
        .method( "clear", &vec::clear )
        // specifically exposing const member functions
        .const_method( "back", &vec::back )
        .const_method( "front", &vec::front )
        .const_method( "at", &vec::at )
        // exposing free functions taking a std::vector < double > * as their first argument
        .method( "assign", &vec_assign )
        .method( "insert", &vec_insert )
        .method( "as.vector", &vec_asR )
        // special methods for indexing
        .const_method( "[[", &vec::at )
        .method( "[[<-", &vec_set )
        ;
}
'

fx_vec <- cxxfunction(signature(), plugin="Rcpp", include=inc_std_vec)
mod_vec <- Module( "mod_vec", getDynLib(fx_vec), mustStart = TRUE )
vec <- mod_vec$vec
v <- new( vec )
v$reserve( 50L )
v$assign( 1:10 )
v$push_back( 10 )
v$size()
v$capacity()
v[[ 0L ]]
v$as.vector()
