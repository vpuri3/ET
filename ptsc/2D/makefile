CFLAGS	         =
FFLAGS	         =
CPPFLAGS         =
FPPFLAGS         =
LOCDIR           = ~/ptsc
EXAMPLESC        = ex1.c
MANSEC           = KSP
CLEANFILES       = 2d 2dtest fftwtest *.o a.out *~ \#*
nNP               = 1

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

2d: 2d.o  chkopts
	-${CLINKER} -o 2d 2d.o  ${PETSC_KSP_LIB} -lfftw3 -lm
	${RM} 2d.o

2dtest: 2dtest.o  chkopts
	-${CLINKER} -o 2dtest 2dtest.o  ${PETSC_KSP_LIB}
	${RM} 2dtest.o

fftwtest: fftwtest.o
	gcc -o fftwtest fftwtest.c -lfftw3 -lm
	${RM} fftwtest.o

#include ${PETSC_DIR}/lib/petsc/conf/test