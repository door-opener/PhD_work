MCC= icc -std=c11 -Wall
#LD= -lgsl -lgslcblas -lm
LDO= -L/home/ucaaitl/gsl/lib/ -lgsl -lgslcblas -lm
LDC= -I/home/uccaitl/gsl/include
OBJ= vanilla_serial.o moveclass.o statmech.o utilities.o statistics.o

vanilla_serial: vanilla_serial.o moveclass.o statmech.o utilities.o statistics.o
	$(MCC) -o vanilla_serial $(OBJ) $(LDO)

vanilla_serial.o: vanilla_serial.c defs.h
	$(MCC) $(LDC) -c vanilla_serial.c

moveclass.o: moveclass.c moveclass.h defs.h
	$(MCC) -c moveclass.c

statistics.o: statistics.c statistics.h defs.h
	$(MCC) $(LDC) -c statistics.c

utilities.o: utilities.c utilities.h defs.h
	$(MCC) -c utilities.c

statmech.o: statmech.c statmech.h defs.h
	$(MCC) $(LDC) -c statmech.c

clean:
	rm -rf *.o
