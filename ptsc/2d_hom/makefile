CFLAGS	         =
FFLAGS	         =
CPPFLAGS         =
FPPFLAGS         =
LOCDIR           = ~/ptsc
EXAMPLESC        = ex1.c
MANSEC           = KSP
CLEANFILES       = 3d *.o a.out *~ \#*
nNP               = 1

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

3d: 3d.o  chkopts
	-${CLINKER} -o 3d 3d.o  ${PETSC_KSP_LIB} -lfftw3 -lm
	${RM} 3d.o


include ${PETSC_DIR}/lib/petsc/conf/test