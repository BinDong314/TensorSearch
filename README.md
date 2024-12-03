# TensorSearch

This repo contains the code developed by the below research paper.
Please drop Bin Dong (dbin@lbl.gov) a message if you need help with the usage of the code.

```
TensorSearch: Parallel Similarity Search on Tensors, IEEE BigData 2024
```

Please see the copyright at the end and License at license.txt

# A simple usage guide
## Install Dependence and TensorSearch
We needs C++, MPI, FasTensor and HDF5 to run TensorSearch.

### C++, MPI
Most systems has C++ installation. Just note that TensorSearch requires "-std=c++17". MPI can be installed from [MPICH](https://www.mpich.org/downloads/) or [OpenMPI](https://docs.open-mpi.org/en/v5.0.x/). In this simple usage guide, I installed MPICH via below command on my MacOS

```shell
brew install mpich
```

### HDF5
Please refer to [HDF5](https://github.com/HDFGroup/hdf5) for details of how to get and install it. Below is my command to compile and install HDF5.

```shell
./configure --prefix=/Users/dbin/work/soft/hdf5-git-dev/build --enable-parallel
```


### FasTensor
Please refer to [FasTensor](https://github.com/BinDong314/FasTensor) for details of how to get and install it. Here, I just compile and install it like this

```shell
./configure --prefix=/Users/dbin/work/fastensor/build-nodash --with-hdf5=/Users/dbin/work/soft/hdf5-git-dev/build CXX=mpic++ CC=mpicc
make 
make install
```

### FasTensor
Change the Makefile to compile the code

```shell
###########################################################
#Please change below configuration for proper installation#
###########################################################
#Compiler of MPI
CCC=mpicxx 
#Installation directory of FasTensor
AU_DIR=/Users/dbin/work/fastensor/build-nodash
#Installation directory of HDF5
HDF5_DIR=/Users/dbin/work/soft/hdf5-git-dev/build
```

After change the AU_DIR and HDF5_DIR to point to right direct, I just call below command to compile the code 
```shell
make
```

## Demo the Usage of TensorSearch
### Generate the data
```shell
./data-generator.sh
```
It produces datasets, including both a DB to search against and a pattern to search with. Now the datasets includes 1D, 2D, and 3D data. Take 2D data as example, it produces three directories and each has a few files. 

```shell
$ ls db-data-2d  # the DB to search against
data-1.h5       data-2.h5       data-3.h5       data-4.h5

$ ls db-pattern-2d # pattern to search with, each file is a pattern
pattern-1.h5    pattern-2.h5    pattern-3.h5    pattern-4.h5

$ ls db-data-pattern-2d # store the result
```
Below is the content of the file db-data-2d/data-1.h5 which is print with HDF5 command h5dump. It is a 2D dataset. The pattern-1.h5 file and files are 2D data too.

```shell
$ h5dump db-data-2d/data-1.h5
HDF5 "db-data-2d/data-1.h5" {
GROUP "/" {
   GROUP "testg" {
      DATASET "testd" {
         DATATYPE  H5T_IEEE_F32LE
         DATASPACE  SIMPLE { ( 4, 4 ) / ( 4, 4 ) }
         DATA {
         (0,0): 7.5, 49.5, 73.5, 58.5,
         (1,0): 30.5, 72.5, 44.5, 78.5,
         (2,0): 23.5, 9.5, 40.5, 65.5,
         (3,0): 92.5, 42.5, 87.5, 3.5
         }}}}}
```

### Run TensorSearch
Below is the command I used to run the search on the 2D dataset above
```shell
./tensor-search -d db-data-2d:/testg/testd  -q db-pattern-2d:/testg/testd -r db-data-pattern-2d:/testg/testd

```

Below is the explanation of the command lines:

```shell
Usage: ./tensor-search [OPTION]
          -h help (--help)
          -d string, [directory/file name]:[dataset name] of the database to search.
          -q string, [directory/file name]:[dataset name]of query/pattern to search with.
          -o string, [directory/file name] to store result. By default similary will be stored in /similarity with /index. 
          -m string, measurement of distance, dotproduct(default), euclidean, cosine, jaccard, angular,  
          -s integer, enable shift with size on database (large tensor) to search, -1 means to shift by the size of pattern 
          -t, enable transport of data from column-wise to row-wise (default input) 
          -k integer, the top k similarity,  result file will be used with [index] dataset to record top-k  
          -r integer, the r-th dimension for the similarity 
          Example: mpirun -n 1 ./tensor-search -d db-data-2d:/testg/testd  -q db-pattern-2d:/testg/testd -r db-data-pattern-2d:/testg/testd
```

****************************

*** Copyright Notice ***

TensorSearch Copyright (c) 2024, The Regents of the University of California,
through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy). All rights reserved.

If you have questions about your rights to use or distribute this software,
please contact Berkeley Lab's Intellectual Property Office at
IPO@lbl.gov.

NOTICE.  This Software was developed under funding from the U.S. Department
of Energy and the U.S. Government consequently retains certain rights.  As
such, the U.S. Government has been granted for itself and others acting on
its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the
Software to reproduce, distribute copies to the public, prepare derivative 
works, and perform publicly and display publicly, and to permit others to do so.


****************************

Please see License Agreement  in license.txt
