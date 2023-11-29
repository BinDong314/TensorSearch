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


export HDF5_USE_FILE_LOCKING=FALSE;$run_command -n 2 ./tensor-search -d db-data-2d/data-1.h5:/testg/testd  -q db-data-2d/data-2.h5:/testg/testd -o test1d-single-file-sm.hdf5 -r 1  -s -1 -p 0 -k 2


# GROUP "/" {
#    DATASET "index" {
#       DATATYPE  H5T_STD_U64LE
#       DATASPACE  SIMPLE { ( 4, 4 ) / ( 4, 4 ) }
#       DATA {
#       (0,0): 0, 0, 0, 0,
#       (1,0): 1, 1, 1, 1,
#       (2,0): 2, 2, 2, 2,
#       (3,0): 3, 3, 3, 3
#       }
#    }
#    DATASET "similarity" {
#       DATATYPE  H5T_IEEE_F32LE
#       DATASPACE  SIMPLE { ( 4, 4 ) / ( 4, 4 ) }
#       DATA {
#       (0,0): 11331, 11680.5, 7455, 9433.5,
#       (1,0): 11680.5, 14329, 8349.5, 10071,
#       (2,0): 7455, 8349.5, 6573, 6350.5,
#       (3,0): 9433.5, 10071, 6350.5, 18031
#       }
#    }
# }
# }