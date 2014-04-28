# For GNU make
FC = ifort
PROJECT = MC3D_Structured_
PROGRAM = Solver
#FFLAGS = -O0 -g -traceback -check all -fp-stack-check -warn all -debug full -module ${OBJ}${PROGRAM}/
FFLAGS = -fast -inline-level=2 -module ${OBJ}${PROGRAM}/
BIN = ../../bin/
OBJ = ../../obj/
INC = -I ../../lib/
LIB = -L ../../lib/

main: RNG.o VAR.o ROUTINES.o IO.o ADV.o main.o
	${FC} ${FFLAGS} -o ${BIN}${PROJECT}${PROGRAM} ${LIB} ${INC} ${OBJ}${PROGRAM}/RNG.o ${OBJ}${PROGRAM}/VAR.o ${OBJ}${PROGRAM}/ROUTINES.o ${OBJ}${PROGRAM}/IO.o ${OBJ}${PROGRAM}/ADV.o ${OBJ}${PROGRAM}/main.o -lMATERIALS_P

%.o: ${PROJECT}%.f90
	${FC} ${FFLAGS} ${LIB} ${INC} -c $< -o ${OBJ}${PROGRAM}/$*.o -lMATERIALS_P