
###########################################################
#Please change below configuration for proper installation#
###########################################################
#Compiler of MPI
CCC=mpicxx 
#Installation directory of FasTensor
AU_DIR=/Users/dbin/work/fastensor/build-nodash
#Installation directory of HDF5
HDF5_DIR=/Users/dbin/work/soft/hdf5-git-dev/build

#Please keep the below unchanged
OTHER_FLAGS=-O3  -std=c++17 
AU_FLAG=-I$(AU_DIR)/include -L$(AU_DIR)/lib -lFastensor
HDF5_FLAG=-I$(HDF5_DIR)/include -L$(HDF5_DIR)/lib -lhdf5  -lhdf5_hl -lz
ALL_FLAGS= $(OTHER_FLAGS) $(AU_FLAG) $(HDF5_FLAG) 

.PHONY:all
all: data-generator tensor-search

data-generator:data-generator.cpp
	$(CCC) -o data-generator data-generator.cpp $(ALL_FLAGS)

tensor-search:tensor-search.cpp
	$(CCC) -o  tensor-search tensor-search.cpp $(ALL_FLAGS) 

clean:
	rm data-generator tensor-search

	
