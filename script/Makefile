CC=g++
CFLAGS=-O3 -std=c++11
LDFLAGS=-static

all: cif2pdb rmligand fasta2nr receptor1 NWalign obo2csv backpropagate mapCSA cif2chain

cif2pdb: cif2pdb.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

rmligand: rmligand.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

fasta2nr: fasta2nr.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

receptor1: receptor1.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

NWalign: NWalign.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

obo2csv: obo2csv.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

backpropagate: backpropagate.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

mapCSA: mapCSA.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

cif2chain: cif2chain.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

