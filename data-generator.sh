#!/bin/bash

DATA_RANK="2"
DATA_SIZE="4,4"
DIR_DB_DATA=db-data-2d
DIR_DB_PATTERN=db-pattern-2d
DIR_DB_DATA_PATTERN=db-data-pattern-2d
mkdir ${DIR_DB_DATA}
mkdir ${DIR_DB_PATTERN}
mkdir ${DIR_DB_DATA_PATTERN}
./data-generator -f ${DIR_DB_PATTERN}/pattern-8x8-1.h5 -g /testg -d /testg/testd -n ${DATA_RANK} -s ${DATA_SIZE} -t 1
./data-generator -f ${DIR_DB_DATA}/data-8x8-1.h5       -g /testg -d /testg/testd -n ${DATA_RANK} -s ${DATA_SIZE} -t 1 -r
cp ${DIR_DB_PATTERN}/pattern-8x8-1.h5  ${DIR_DB_PATTERN}/pattern-8x8-2.h5
cp ${DIR_DB_PATTERN}/pattern-8x8-1.h5  ${DIR_DB_PATTERN}/pattern-8x8-3.h5
cp ${DIR_DB_PATTERN}/pattern-8x8-1.h5  ${DIR_DB_PATTERN}/pattern-8x8-4.h5
cp ${DIR_DB_DATA}/data-8x8-1.h5   ${DIR_DB_DATA}/data-8x8-2.h5
cp ${DIR_DB_DATA}/data-8x8-1.h5   ${DIR_DB_DATA}/data-8x8-3.h5
cp ${DIR_DB_DATA}/data-8x8-1.h5   ${DIR_DB_DATA}/data-8x8-4.h5