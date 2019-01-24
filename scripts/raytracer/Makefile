f77=gfortran

.f.o:
	$(f77) -c $*.f -ffixed-line-length-80 --std=legacy

#  path of executables, leave empty for current dir
BIN =

RAYP=rayP

all: rayP

rayP: rayP.o
	$(f77) -o $(BIN)$(RAYP) $(RAYP).o

clean:
	rm -f *.o *.so rayP rayS
