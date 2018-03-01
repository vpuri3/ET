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

2d: 2d.o  chkopts
	-${CLINKER} -o 2d 2d.o  ${PETSC_KSP_LIB}
	${RM} 2d.o

2dtest: 2dtest.o  chkopts
	-${CLINKER} -o 2dtest 2dtest.o  ${PETSC_KSP_LIB}
	${RM} 2dtest.o


#include ${PETSC_DIR}/lib/petsc/conf/test