CXXFLAGS=-DNDEBUG -std=c++11 -O3 -w -fPIC
LPFLAGS=-Dlpv
CC=g++
#CC=x86_64-w64-mingw32-g++
OBJ=src/ParseCommandLine.o src/common_utils.o 
OBJ2=CSSR/cssr_struct.o CSSR/PDBParser.o 
LDFLAGS=-lstdc++ #-static

exe/CSSR: ${OBJ} CSSR/CSSR.o ${OBJ2}
	${CC} -o exe/CSSR ${CXXFLAGS} CSSR/CSSR.o ${OBJ} ${OBJ2} ${LDFLAGS}

CSSR/CSSR.o: CSSR/CSSR.cpp
	${CC} -c -o $@ ${CXXFLAGS} CSSR/CSSR.cpp
CSSR/PDBParser.o: CSSR/PDBParser.cpp
	${CC} -c -o $@ ${CXXFLAGS} CSSR/PDBParser.cpp
CSSR/cssr_struct.o: CSSR/cssr_struct.cpp
	${CC} -c -o $@ ${CXXFLAGS} CSSR/cssr_struct.cpp ${LPFLAGS}
src/ParseCommandLine.o: src/ParseCommandLine.cpp
	${CC} -c -o $@ ${CXXFLAGS} src/ParseCommandLine.cpp
src/common_utils.o: src/common_utils.cpp
	${CC} -c -o $@ ${CXXFLAGS} src/common_utils.cpp
install: exe/CSSR 
	cp exe/CSSR ../script/CSSR
clean:
	rm ${OBJ} ${OBJ2} exe/CSSR CSSR/CSSR.o 
