CFLAGS	         =
FFLAGS	         =
CPPFLAGS         =
FPPFLAGS         =
LOCDIR           = ~/ptsc
EXAMPLESC        = ex1.c
MANSEC           = KSP
CLEANFILES       = rhs.vtk solution.vtk test *~ \#*\#
nNP               = 1

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

ex1: ex1.o  chkopts
	-${CLINKER} -o ex1 ex1.o  ${PETSC_KSP_LIB}
	${RM} ex1.o

test: test.o  chkopts
	-${CLINKER} -o test test.o  ${PETSC_KSP_LIB}
	${RM} test.o

#include ${PETSC_DIR}/lib/petsc/conf/test