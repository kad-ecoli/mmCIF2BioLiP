CC=g++
CFLAGS=-O3 -std=c++11
LDFLAGS=-static

all: cif2pdb rmligand

cif2pdb: cif2pdb.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

rmligand: rmligand.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}
