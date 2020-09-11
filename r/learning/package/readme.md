## Package structure
file:
* DESCRIPTION
* NAMESPACE
* INDEX
* configure - shell scripts run before install
* cleanup - shell scripts run after install
* LICENSE

subdirectories:
* R - R code files
* data
* demo
* exec - contain additional executable scripts the package needs
* inst - copied recursively to the installation directory
* man
* po -  files related to localization
* src
* tests - additional package-specific test code
* tools - auxiliary files needed during configuration
* vignettes

### build
```sh
R CMD command --help  # command can be (BATCH|COMPILE|SHLIB|INSTALL|build|check|config|...)
pkg-config -h | less   # Return metainformation about installed libraries
R CMD config --help | less
R CMD INSTALL --clean pkg           # make clean
```

tools package
```R
```

### dyn.load
* done automatically based on the useDynLib() declaration in the NAMESPACE file
* On loading a so, R will look for a routine within that so named R_init_lib where lib is the name of the so file with the extension removed

```c++
// src/init.c
#include <R_ext/Rdynload.h>

void myC(double *x, int *n, char **names, int *status);
static R_NativePrimitiveArgType myC_t[] = {
    REALSXP, INTSXP, STRSXP, LGLSXP
};
static const R_CMethodDef cMethods[] = {
   {"myC", (DL_FUNC) &myC, 4, myC_t},
   {NULL, NULL, 0, NULL}
};

SEXP myCall(SEXP a, SEXP b, SEXP c);

static const R_CallMethodDef callMethods[]  = {
  {"myCall", (DL_FUNC) &myCall, 3},
  {NULL, NULL, 0}
};

void R_init_mylib(DllInfo *info) {
  /* Register routines, allocate resources. */
  R_registerRoutines(dll, cMethods/*R_CMethodDef*/, callMethods/*R_CallMethodDef*/, NULL/*R_FortranMethodDef*/, NULL/*R_ExternalMethodDef*/);
  R_useDynamicSymbols(dll, FALSE); // only find registered symbols, others will not to be searched
  
  /**
   * useDynLib(mypkg, .registration = TRUE, .fixes = "C_")
   * now only .Call(C_reg) will work (and not .Call("reg")) 
   */
  R_forceSymbols(dll, TRUE);  
}

void R_unload_mylib(DllInfo *info) {
  /* Release resources. */
}
```

extract cpp declarations
```sh
cproto -I/usr/share/R/include -e student_export.cpp
```

##### link
PKG_LIBS = -L$(XML_DIR)/lib -lxml2 -pthread
PKG_CPPFLAGS = -pthread

### src
sources and headers, and optionally a file Makevars or Makefile

### man subdirectory
documentation files Rd format.
* macros subdirectory contain source for user-defined Rd macros

### NAMESPACE
* specify which variables should be exported to make them available 
* which variables should be imported from other packages
* only the exported variables are placed in the attached frame
* imports variables from other packages will cause these other packages to be loaded as well

```R
import(pkg1, pkg2)
importFrom(pkg1, var1, var2)
import(pkg1, except=c(var3, var4))  # imports every symbol from pkg1 except var3 and var4
```

#### useDynLib
allows shared objects that need to be loaded

* useDynLib(foo, myRoutine) enables .Call(myRoutine, x, y) instead of .Call("myRoutine", x, y, PACKAGE = "foo")
* useDynLib(myDLL, .registration = TRUE) causes the so to be loaded and also for the R variables of all which has registration information to be defined in the package’s namespace

#### S4 class export/inport
* exportClasses / exportMethods / exportClassPattern
* importClassesFrom(package, ...) / importMethodsFrom(package, ...)

### DESCRIPTION
mandatory fields = 'Package', 'Version', 'License', 'Description', 'Title', 'Author', 'Maintainer'
* 'Depends' gives a list of package names which this package depends on
* 'Imports' lists packages whose namespaces are imported from (as specified in the NAMESPACE file) but which do not need to be attached
* 'LinkingTo' make use of header files in other packages to compile its cpp code

### data subdirectory
three types 
* plain R code (.R or .r)
* tables (.tab, .txt, or .csv, see ?data for the file formats, and note that .csv is not the standard CSV format)
* save() images (.RData or .rda)

### demo subdirectory
R scripts (for running via demo()) that demonstrate some of the functionality of the package

### compile C and call
```c++
// foo.c
#include <Rinternals.h>

SEXP foo_func(SEXP x) {
    return x;
}
```

compile it:
```bash
R CMD SHLIB foo.c   # compile to foo.so
R CMD SHLIB --help
# PKG_CPPFLAGS  - mainly ‘-I’, ‘-D’ and ‘-U’ flags
# PKG_LIBS      - additional ‘-l’ and ‘-L’ flags to linker
```

wrapper it in R:
```R
# wrappers.R 
# wrapper function to invoke helloA1
dyn.load("foo.so")
foo2 <- function(x) .Call("foo_func", x, PACKAGE = "foo")
```

use it in client code:
```R
source('wrappers.R')
greeting <- helloA1()
```

## Document
### RD format
src/library/${pkg_name}/man/${func_name}.Rd
* function sections, (name|alias|title|description|usage|arguments|seealso|examples|keyword)

process:
```sh
R CMD Rd2pdf
R CMD Sweave
R CMD Stangle
R CMD Rdconv --help
```

function example:
```func Rd
% comment, function load document, src/library/base/man/load.Rd
\name{load}
\alias{load}
\title{func title}
\description{func desc}
\usage{
func(arg1, arg2, …)
}
\arguments{
  \item{arg_i}{Description of arg_i}
}
\seealso{
  \code{\link{related_func}}.
}
\examples{
## save all data
save(list = ls(), file= "all.RData")
}
\keyword{file}
```

data set example:
```dataset Rd
% src/library/datasets/man/rivers.Rd
\name{rivers}
\docType{data}
\alias{rivers}
\title{dataset title}
\description{data set desc}
\usage{rivers}
\format{A vector containing 141 observations.}
\source{World Almanac and Book of Facts, 1975, page 406.}
\keyword{datasets}
```

S4 classes and methods
```class Rd
\S4method{generic}{signature_list}(argument_list)
\S4method{coerce}{ANY,NULL}(from, to)
```