apt install libopenblas-dev liblapack-dev libarpack2-dev libsuperlu-dev
cd armadillo-9.900.1/
./configure
# disable hdf5
cmake -D DETECT_HDF5=false .
make
sudo make install

# compile with mkl
cpu-info                                # check if intel cpu
sh l_onemkl_p_2021.1.1.52_offline.sh
ln -s /opt/intel/oneapi/mkl/2021.1.1 /opt/intel/mkl
vim /etc/ld.so.conf.d/libmkl.conf  # add /opt/intel/mkl/lib/intel64/
ldconfig