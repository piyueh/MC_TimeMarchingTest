Test: RNG.o MC_TimeMarchingTest_Variables.o MC_TimeMarchingTest_Adv.o MC_TimeMarchingTest_Main.o
	ifort -traceback -o ./bin/MC_TimeMarchingTest ./obj/RNG.o ./obj/MC_TimeMarchingTest_Variables.o ./obj/MC_TimeMarchingTest_Adv.o ./obj/MC_TimeMarchingTest_Main.o

RNG.o: ./SourceCodes/RNG.f90
	ifort -traceback -warn stderrors -warn errors -warn all -debug full -module ./obj/ -c ./SourceCodes/RNG.f90 -o ./obj/RNG.o

MC_TimeMarchingTest_Variables.o: ./SourceCodes/MC_TimeMarchingTest_Variables.f90 RNG.o
	ifort -traceback -warn stderrors -warn errors -warn all -debug full -module ./obj/ -c ./SourceCodes/MC_TimeMarchingTest_Variables.f90 -o ./obj/MC_TimeMarchingTest_Variables.o

MC_TimeMarchingTest_Adv.o: ./SourceCodes/MC_TimeMarchingTest_Adv.f90 RNG.o MC_TimeMarchingTest_Variables.o
	ifort -traceback -warn stderrors -warn errors -warn all -debug full -module ./obj/ -c ./SourceCodes/MC_TimeMarchingTest_Adv.f90 -o ./obj/MC_TimeMarchingTest_Adv.o

MC_TimeMarchingTest_Main.o: ./SourceCodes/MC_TimeMarchingTest_Main.f90 RNG.o MC_TimeMarchingTest_Variables.o MC_TimeMarchingTest_Adv.o
	ifort -traceback -warn stderrors -warn errors -warn all -debug full -module ./obj/ -c ./SourceCodes/MC_TimeMarchingTest_Main.f90 -o ./obj/MC_TimeMarchingTest_Main.o

clean:
	rm ./obj/*
	rm ./bin/*