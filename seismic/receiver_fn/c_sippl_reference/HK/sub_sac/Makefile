CFLAGS = -g -c -Wall -O2 -funroll-loops -fomit-frame-pointer # -DSOLARIS -DDEBUG

CC = g++
FF = gfortran

LIBOBJ = sac.o	
LIB    = ../lib/libsaccp.a


all: $(LIBOBJ) arx lib

arx:
	ar -q $(LIB)

lib: $(LIBOBJ)
	ar ruv $(LIB) $(LIBOBJ)
	ranlib $(LIB)

sac.o: sac.cpp sac.h
	$(CC) $(CFLAGS) sac.cpp 

.PHONY: clean
clean:
	\rm *.o
