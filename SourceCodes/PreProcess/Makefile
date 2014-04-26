# For GNU make
FC = ifort
PROJECT = MC3D_Structured
PROGRAM = PreProcessing/
#FFLAGS = -O0 -g -traceback -check all -fp-stack-check -warn all -debug full -module ${OBJ}${PROGRAM}
FFLAGS = -fast -inline-level=2 -module ${OBJ}${PROGRAM}
BIN = ../../bin/
OBJ = ../../obj/
INC = -I ../../lib/
LIB = -L ../../lib/

all:
	if [ ! -e ${BIN} ]; then mkdir ${BIN}; fi
	if [ ! -e ${OBJ} ]; then mkdir ${OBJ}; fi
	if [ ! -e ${BIN}${PROGRAM} ]; then mkdir ${BIN}${PROGRAM}; fi
	if [ ! -e ${OBJ}${PROGRAM} ]; then mkdir ${OBJ}${PROGRAM}; fi
	make clean
	make main

main: VAR.o PreProcessing.o
	${FC} ${FFLAGS} -o ${BIN}${PROGRAM}${PROJECT} ${LIB} ${INC} ${OBJ}${PROGRAM}VAR.o ${OBJ}${PROGRAM}PreProcessing.o -lMATERIALS

%.o: ${PROJECT}_%.f90
	${FC} ${FFLAGS} ${LIB} ${INC} -c $< -o ${OBJ}${PROGRAM}$*.o -lMATERIALS

.PHONY: clean
clean:
	@ - rm ${OBJ}${PROGRAM}*
	@ - rm ${BIN}${PROGRAM}*