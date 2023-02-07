CC=g++
MINGW=x86_64-w64-mingw32-g++ -static
CFLAGS=-O3 -ffast-math
LDFLAGS=#-static# -lm
PROGRAM=USalign pdb2xyz pdb2fasta xyz_sfetch

all: ${PROGRAM}

USalign: USalign.cpp SOIalign.h MMalign.h param_set.h basic_fun.h Kabsch.h NW.h TMalign.h pstream.h se.h NWalign.h BLOSUM.h flexalign.h
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

pdb2xyz: pdb2xyz.cpp basic_fun.h pstream.h
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

xyz_sfetch: xyz_sfetch.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

pdb2fasta: pdb2fasta.cpp basic_fun.h pstream.h
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

clean:
	rm -f ${PROGRAM}

install: all
	cp ${PROGRAM} ../script/
