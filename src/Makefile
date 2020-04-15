

FC=ftn -extend-source
OPT=-g -check all -fpe0
EXE=../iotest

OBJ=  main.o setups.o bindmp.o
${EXE}: ${OBJ}
	${FC} ${OPT} main.o bindmp.o setups.o -o ${EXE}

clean:
	rm -fr *~ *.o *.mod ${EXE}

.SUFFIXES: .F.o

.F.o:
	${FC} ${OPT} -c $*.F

main.o: setups.o bindmp.o

setups.o:

bindmp.o: setups.o

