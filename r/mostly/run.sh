Rscript script.r
Rscript -e 'install.packages(c("devtools", "Rcpp"))'

R CMD config --cppflags
R CMD build .
R CMD build pkg
R CMD INSTALL pkg_1.0.tar.gz
R CMD INSTALL --build pkg
R CMD INSTALL --clean pkg           # make clean
R CMD check
R CMD check --help

# debug script
R -d gdb --vanilla
R -d valgrind --vanilla < mypkg.R
R -d "valgrind --tool=memcheck --leak-check=full" --vanilla < mypkg.R