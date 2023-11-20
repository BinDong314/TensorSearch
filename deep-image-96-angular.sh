#!/bin/bash

# Initialize the MPI command variable
run_command=""

# Check if mpirun is available
if command -v mpirun &> /dev/null; then
    run_command="mpirun"
elif command -v srun &> /dev/null; then
    run_command="srun"
else
    echo "Error: Neither mpirun nor srun is available."
    exit 1
fi


export HDF5_USE_FILE_LOCKING=FALSE;$run_command -n 2 ./tensor-search -d ./compare/deep-image-96-angular.hdf5:/train  -q ./compare/deep-image-96-angular.hdf5:/test -r ./compare/deep-image-96-angular-ts-result.hdf5:/similarity -k 1  -s -1 -p 0
