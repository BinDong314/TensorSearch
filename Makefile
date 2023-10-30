#Bin's Mac

AU_DIR=/Users/dbin/work/fastensor/build
HDF5_DIR=/Users/dbin/work/soft/hdf5-git-dev/build
DASH_DIR=/Users/dbin/work/soft/dash/build/install

CCC=mpicxx
OTHER_FLAGS=-O3  -std=c++17 -I./DasLib/

AU_FLAG=-I$(AU_DIR)/include -L$(AU_DIR)/lib -lFastensor
HDF5_FLAG=-I$(HDF5_DIR)/include -L$(HDF5_DIR)/lib -lhdf5  -lhdf5_hl -lz
DASH_FLAG=-I$(DASH_DIR)/include -L$(DASH_DIR)/lib -ldash-mpi -ldart-mpi -ldart-base -lpthread -DDASH_ENABLE_HDF5 -DDASH_MPI_IMPL_ID='mpich' -DHAS_DASH_ENDPOINT
ALL_FLAGS= $(OTHER_FLAGS) $(AU_FLAG) $(HDF5_FLAG) $(DASH_FLAG)
#$(FFTW_FLAG)  $(DASH_FLAG) $(EIGEN3_FLAG) 

.PHONY:all
all: data-generator tensor-search

data-generator:data-generator.cpp
	$(CCC) -o data-generator data-generator.cpp $(ALL_FLAGS)

tensor-search:tensor-search.cpp
	$(CCC) -o  tensor-search tensor-search.cpp $(ALL_FLAGS) 

clean:
	rm data-generator tensor-search

	
