CC=g++
CFLAGS=-O3 -std=c++11
LDFLAGS=-static

all: cif2pdb rmligand fasta2nr

cif2pdb: cif2pdb.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

rmligand: rmligand.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

fasta2nr: fasta2nr.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}
