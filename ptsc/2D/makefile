CFLAGS	         =
FFLAGS	         =
CPPFLAGS         =
FPPFLAGS         =
LOCDIR           = ~/ptsc
EXAMPLESC        = ex1.c
MANSEC           = KSP
CLEANFILES       = 2d test fftw1d fftw2d *.o a.out *~ \#*
nNP               = 1

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

# hey Pablo you should be able to see this.

2d: 2d.o  chkopts
	-${CLINKER} -o 2d 2d.o  ${PETSC_KSP_LIB} -lfftw3 -lm
	${RM} 2d.o

test: test.o  chkopts
	-${CLINKER} -o test test.o  ${PETSC_KSP_LIB}
	${RM} test.o

fftw1d: fftw1d.o
	gcc -o fftw1d fftw1d.c -lfftw3 -lm
	${RM} fftw1d.o

fftw2d: fftw2d.o
	gcc -o fftw2d fftw2d.c -lfftw3 -lm
	${RM} fftw2d.o

#include ${PETSC_DIR}/lib/petsc/conf/test