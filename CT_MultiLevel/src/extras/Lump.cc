plist = (CCTK_REAL**) realloc(plist, 3*sizeof(CCTK_REAL*));

plist[0] = (CCTK_REAL*) CCTK_VarDataPtr(cctkGH, 0, "CT_Analytic::testdxK");
plist[1] = (CCTK_REAL*) CCTK_VarDataPtr(cctkGH, 0, "CT_Analytic::testdyK");
plist[2] = (CCTK_REAL*) CCTK_VarDataPtr(cctkGH, 0, "CT_Analytic::testdzK");
