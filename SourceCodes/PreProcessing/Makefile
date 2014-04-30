# For GNU make
FC = ifort
PROJECT = MC3D_Structured_
PROGRAM = PreProcessing
BIN = ../../bin/
OBJ = ../../obj/
INC = -I ../../lib/
LIB = -L ../../lib/


debug: FFLAGS = -O0 -g -traceback -check all -fp-stack-check -warn all -debug full -module ${OBJ}${PROGRAM}/
debug: VAR.o PreProcessing.o
	${FC} ${FFLAGS} -o ${BIN}${PROJECT}${PROGRAM} ${LIB} ${OBJ}${PROGRAM}/VAR.o ${OBJ}${PROGRAM}/PreProcessing.o -lMATERIALS

release: FFLAGS = -fast -inline-level=2 -module ${OBJ}${PROGRAM}/
release: VAR.o PreProcessing.o
	${FC} ${FFLAGS} -o ${BIN}${PROJECT}${PROGRAM} ${LIB} ${OBJ}${PROGRAM}/VAR.o ${OBJ}${PROGRAM}/PreProcessing.o -lMATERIALS

%.o: ${PROJECT}%.F90
	${FC} ${FFLAGS} ${LIB} ${INC} -c $< -o ${OBJ}${PROGRAM}/$*.o -lMATERIALS

clean:
	- @ rm ${OBJ}${PROGRAM}/*
	- @ rm ${BIN}${PROJECT}${PROGRAM}