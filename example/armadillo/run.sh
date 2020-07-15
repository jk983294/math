apt install libopenblas-dev liblapack-dev libarpack2-dev libsuperlu-dev
cd armadillo-9.900.1/
./configure
# disable hdf5
cmake -D DETECT_HDF5=false .
make
sudo make install