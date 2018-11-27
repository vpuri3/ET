CFLAGS	         =
FFLAGS	         =
CPPFLAGS         =
FPPFLAGS         =
LOCDIR           = ~/ptsc
EXAMPLESC        = ex1.c
MANSEC           = KSP
CLEANFILES       = 2d *.o a.out *~ \#*
nNP               = 1

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

2d: 2d.o  chkopts
	-${CLINKER} -o 2d 2d.o  ${PETSC_KSP_LIB} -lfftw3 -lm
	${RM} 2d.o


include ${PETSC_DIR}/lib/petsc/conf/test