CFLAG = -g -c -Wall -O2 -funroll-loops -fomit-frame-pointer # -DSOLARIS -DDEBUG

CC = g++
FF = gfortran

ALLOBJ = util_sub.o tseris_proc.o util_time.o distaz.o
LIB    = ../lib/libutil.a
INCSAC = -I../sub_sac
INCNR = -I../sub_nrc


all: $(ALLOBJ) arx lib

arx:
	ar -q $(LIB)

lib: $(ALLOBJ)
	ar ruv $(LIB) $(ALLOBJ)
	ranlib $(LIB)

util_sub.o: util_sub.cpp
	$(CC) $(CFLAG) util_sub.cpp

tseris_proc.o: tseris_proc.cpp
	$(CC) $(CFLAG) tseris_proc.cpp $(INCNR)

util_time.o: util_time.cpp util_time.h
	$(CC) $(CFLAG) util_time.cpp

distaz.o: distaz.cpp
	$(CC) $(CFLAG) distaz.cpp


.PHONY: clean
clean:
	\rm *.o
