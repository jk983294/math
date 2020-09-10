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
```c
// helloA1.c
#include <R.h>
 
SEXP helloA1() {
  printf("Hello World!\n");
  return(R_NilValue);
}
```

```bash
R CMD SHLIB helloA1.c   # compile to helloA1.so
```

```R
# wrappers.R 
# wrapper function to invoke helloA1
dyn.load("helloA1.so")
helloA1 <- function() {
  result <- .Call("helloA1")
}
```

```R
source('wrappers.R')
greeting <- helloA1()
```