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


export HDF5_USE_FILE_LOCKING=FALSE;$run_command -n 1 ./tensor-search -d app/absj-sub.h5:/absj -q app/absj-sub-template.h5:/absj -o app/absj-sub-results.h5  -s -1 
