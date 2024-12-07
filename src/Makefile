##########################
# setup for XC50 of CfCA
##########################
FC=mpiifx -fc=ifort -diag-disable=10448
#OPT=-g -O0 -check all -fpe0
#LOPT=-Wl,--no-relax
OPT=-g -O2

##########################
# setup for Msize Server
##########################
#FC=mpif90
#OPT=-g

##########################
# Main
##########################

EXE=../iotest

OBJ=  main.o setups.o bindmp.o

${EXE}: ${OBJ}
	${FC} ${OPT} ${LOPT}  main.o bindmp.o setups.o -o ${EXE}

.SUFFIXES: .f90 .o

.f90.o:
	${FC} ${OPT} -c $*.f90

main.o: setups.o bindmp.o

setups.o:

bindmp.o: setups.o

##########################
# PHONY
##########################

.PHONY: clean

clean:
	rm -fr *~ *.o *.mod ${EXE}

