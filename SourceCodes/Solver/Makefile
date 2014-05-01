# For GNU make
FC = ifort
PROJECT = MC3D_Structured_
PROGRAM = Solver
BIN = ../../bin/
OBJ = ../../obj/
INC = -I ../../lib/
LIB = -L ../../lib/

debug: FFLAGS = -O0 -g -traceback -check all -fp-stack-check -warn all -debug full -module ${OBJ}${PROGRAM}/
debug: RNG.o VAR.o ROUTINES.o ADV.o HEAT.o IO.o main.o
	${FC} ${FFLAGS} -o ${BIN}${PROJECT}${PROGRAM} ${LIB} ${OBJ}${PROGRAM}/RNG.o ${OBJ}${PROGRAM}/VAR.o ${OBJ}${PROGRAM}/ROUTINES.o ${OBJ}${PROGRAM}/IO.o ${OBJ}${PROGRAM}/ADV.o ${OBJ}${PROGRAM}/HEAT.o ${OBJ}${PROGRAM}/main.o -lMATERIALS

release: FFLAGS = -fast -inline-level=2 -module ${OBJ}${PROGRAM}/
release: RNG.o VAR.o ROUTINES.o ADV.o HEAT.o IO.o main.o
	${FC} ${FFLAGS} -o ${BIN}${PROJECT}${PROGRAM} ${LIB} ${OBJ}${PROGRAM}/RNG.o ${OBJ}${PROGRAM}/VAR.o ${OBJ}${PROGRAM}/ROUTINES.o ${OBJ}${PROGRAM}/IO.o ${OBJ}${PROGRAM}/ADV.o ${OBJ}${PROGRAM}/HEAT.o ${OBJ}${PROGRAM}/main.o -lMATERIALS

%.o: ${PROJECT}%.F90
	${FC} ${FFLAGS} ${LIB} -c $< -o ${OBJ}${PROGRAM}/$*.o -lMATERIALS

.PHONY: clean
clean:
	- @ rm ${OBJ}${PROGRAM}/*
	- @ rm ${BIN}${PROJECT}${PROGRAM}