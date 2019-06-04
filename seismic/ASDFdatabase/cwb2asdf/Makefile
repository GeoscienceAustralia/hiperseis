all: demultiplex

CC=g++
CFLAGS=-O3
LIB_PATH=../../../libmseed
INC_PATH=../../../libmseed

SRCS=demultiplex.cc

clean:
	rm -f demultiplex

.PHONY:demultiplex

demultiplex:
	${CC} -o demultiplex ${SRCS} ${CFLAGS} -L${LIB_PATH} -I${INC_PATH} -lmseed
