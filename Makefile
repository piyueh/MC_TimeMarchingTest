# For GNU make
FC = ifort
PROJECT = MC_TimeMarchingTest
#FFLAGS = -O0 -g -traceback -check all -fp-stack-check -warn all -debug full -module ${OBJ}
FFLAGS = -O3 -ipo -fast -inline-level=2 -module ${OBJ}
SRC = ./SourceCodes/
BIN = ./bin/
OBJ = ./obj/
INC = -I
LIB = -L

all:
	if [ ! -e ${BIN} ]; then mkdir ${BIN}; fi
	if [ ! -e ${OBJ} ]; then mkdir ${OBJ}; fi
	make clean
	make main

main: RNG.o VAR.o ROUTINES.o INIT.o ADV.o main.o
	${FC} ${FFLAGS} -o ${BIN}${PROJECT} ./obj/RNG.o ./obj/VAR.o ./obj/ROUTINES.o ./obj/INIT.o ./obj/ADV.o ./obj/main.o

%.o: ${SRC}${PROJECT}_%.f90
	${FC} ${FFLAGS} -c $< -o ${OBJ}$*.o

.PHONY: clean
clean:
	@ - rm ${OBJ}*
	@ - rm ${BIN}*