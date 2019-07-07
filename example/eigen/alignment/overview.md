### dixed-size vectorizable Eigen objects
An Eigen object is called "fixed-size vectorizable" if it has fixed size and that size is a multiple of 16 bytes.

Eigen will request 16-byte alignment for these objects, and henceforth rely on these objects being aligned so no runtime check for alignment is performed.

### Passing Eigen objects by value to functions
passing fixed-size vectorizable Eigen objects by value is not only inefficient, it can be illegal or make your program crash! And the reason is that these Eigen objects have alignment modifiers that aren't respected when they are passed by value.

please pass by reference

there is no problem with functions that return objects by value.
