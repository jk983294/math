# https://colinfay.me/writing-r-extensions/debugging.html

R -d valgrind --vanilla < test.R
R -d "valgrind --tool=memcheck --leak-check=full" --vanilla < mypkg-Ex.R