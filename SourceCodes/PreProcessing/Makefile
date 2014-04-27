# For GNU make
FC = ifort
PROJECT = MC3D_Structured_
PROGRAM = PreProcessing
#FFLAGS = -O0 -g -traceback -check all -fp-stack-check -warn all -debug full -module ${OBJ}${PROGRAM}/
FFLAGS = -fast -inline-level=2 -module ${OBJ}${PROGRAM}/
BIN = ../../bin/
OBJ = ../../obj/
INC = -I ../../lib/
LIB = -L ../../lib/

main: VAR.o PreProcessing.o
	${FC} ${FFLAGS} -o ${BIN}${PROJECT}${PROGRAM} ${LIB} ${INC} ${OBJ}${PROGRAM}/VAR.o ${OBJ}${PROGRAM}/PreProcessing.o -lMATERIALS

%.o: ${PROJECT}%.f90
	${FC} ${FFLAGS} ${LIB} ${INC} -c $< -o ${OBJ}${PROGRAM}/$*.o -lMATERIALS