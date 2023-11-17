#!/bin/bash

DATA_RANK="2"
DATA_SIZE="4,4"
DIR_DB_DATA=db-data-2d
DIR_DB_PATTERN=db-pattern-2d
DIR_DB_DATA_PATTERN=db-data-pattern-2d
mkdir ${DIR_DB_DATA}
rm ${DIR_DB_DATA}/*
mkdir ${DIR_DB_PATTERN}
rm ${DIR_DB_PATTERN}/*
mkdir ${DIR_DB_DATA_PATTERN}
rm ${DIR_DB_DATA_PATTERN}/*
./data-generator -f ${DIR_DB_PATTERN}/pattern-1.h5 -g /testg -d /testg/testd -n ${DATA_RANK} -s ${DATA_SIZE} -t 1
./data-generator -f ${DIR_DB_DATA}/data-1.h5       -g /testg -d /testg/testd -n ${DATA_RANK} -s ${DATA_SIZE} -t 1 -r
cp ${DIR_DB_PATTERN}/pattern-1.h5  ${DIR_DB_PATTERN}/pattern-2.h5
cp ${DIR_DB_PATTERN}/pattern-1.h5  ${DIR_DB_PATTERN}/pattern-3.h5
cp ${DIR_DB_PATTERN}/pattern-1.h5  ${DIR_DB_PATTERN}/pattern-4.h5
cp ${DIR_DB_DATA}/data-1.h5   ${DIR_DB_DATA}/data-2.h5
cp ${DIR_DB_DATA}/data-1.h5   ${DIR_DB_DATA}/data-3.h5
cp ${DIR_DB_DATA}/data-1.h5   ${DIR_DB_DATA}/data-4.h5




DATA_RANK="1"
DATA_SIZE="4"
DIR_DB_DATA=db-data-1d
DIR_DB_PATTERN=db-pattern-1d
DIR_DB_DATA_PATTERN=db-data-pattern-1d
mkdir ${DIR_DB_DATA}
rm ${DIR_DB_DATA}/*
mkdir ${DIR_DB_PATTERN}
rm ${DIR_DB_PATTERN}/*
mkdir ${DIR_DB_DATA_PATTERN}
rm ${DIR_DB_DATA_PATTERN}/*
./data-generator -f ${DIR_DB_PATTERN}/pattern-1.h5 -g /testg -d /testg/testd -n ${DATA_RANK} -s ${DATA_SIZE} -t 1
./data-generator -f ${DIR_DB_DATA}/data-1.h5       -g /testg -d /testg/testd -n ${DATA_RANK} -s ${DATA_SIZE} -t 1 -r
cp ${DIR_DB_PATTERN}/pattern-1.h5  ${DIR_DB_PATTERN}/pattern-2.h5
cp ${DIR_DB_PATTERN}/pattern-1.h5  ${DIR_DB_PATTERN}/pattern-3.h5
cp ${DIR_DB_PATTERN}/pattern-1.h5  ${DIR_DB_PATTERN}/pattern-4.h5
cp ${DIR_DB_DATA}/data-1.h5   ${DIR_DB_DATA}/data-2.h5
cp ${DIR_DB_DATA}/data-1.h5   ${DIR_DB_DATA}/data-3.h5
cp ${DIR_DB_DATA}/data-1.h5   ${DIR_DB_DATA}/data-4.h5




DATA_RANK="3"
DATA_SIZE="4,4,4"
DIR_DB_DATA=db-data-3d
DIR_DB_PATTERN=db-pattern-3d
DIR_DB_DATA_PATTERN=db-data-pattern-3d
mkdir ${DIR_DB_DATA}
rm ${DIR_DB_DATA}/*
mkdir ${DIR_DB_PATTERN}
rm ${DIR_DB_PATTERN}/*
mkdir ${DIR_DB_DATA_PATTERN}
rm ${DIR_DB_DATA_PATTERN}/*
./data-generator -f ${DIR_DB_PATTERN}/pattern-1.h5 -g /testg -d /testg/testd -n ${DATA_RANK} -s ${DATA_SIZE} -t 1
./data-generator -f ${DIR_DB_DATA}/data-1.h5       -g /testg -d /testg/testd -n ${DATA_RANK} -s ${DATA_SIZE} -t 1 -r
cp ${DIR_DB_PATTERN}/pattern-1.h5  ${DIR_DB_PATTERN}/pattern-2.h5
cp ${DIR_DB_PATTERN}/pattern-1.h5  ${DIR_DB_PATTERN}/pattern-3.h5
cp ${DIR_DB_PATTERN}/pattern-1.h5  ${DIR_DB_PATTERN}/pattern-4.h5
cp ${DIR_DB_DATA}/data-1.h5   ${DIR_DB_DATA}/data-2.h5
cp ${DIR_DB_DATA}/data-1.h5   ${DIR_DB_DATA}/data-3.h5
cp ${DIR_DB_DATA}/data-1.h5   ${DIR_DB_DATA}/data-4.h5