# For GNU make
SRC = ./SourceCodes/
PROGRAM1 = PreProcessing
PROGRAM2 = Solver
BIN = ./bin/
OBJ = ./obj/

all:
	make clean
	if [ ! -e ${BIN} ]; then mkdir ${BIN}; fi
	if [ ! -e ${OBJ} ]; then mkdir ${OBJ}; fi
	if [ ! -e ${OBJ}${PROGRAM1}/ ]; then mkdir ${OBJ}${PROGRAM1}/; fi
	if [ ! -e ${OBJ}${PROGRAM2}/ ]; then mkdir ${OBJ}${PROGRAM2}/; fi
	make PreProcessing
	make Solver

PreProcessing:
	cd ${SRC}${PROGRAM1}/ && make main

Solver:
	cd ${SRC}${PROGRAM2}/ && make main


.PHONY: clean
clean:
	@ - rm -r ${OBJ}/*
	@ - rm -r ${BIN}/*