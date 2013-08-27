plist = (CCTK_REAL**) realloc(plist, 5*sizeof(CCTK_REAL*));

plist[0] = (CCTK_REAL*) CCTK_VarDataPtr(cctkGH, 0, "CT_Analytic::testdxK");
plist[1] = (CCTK_REAL*) CCTK_VarDataPtr(cctkGH, 0, "CT_Analytic::testdyK");
plist[2] = (CCTK_REAL*) CCTK_VarDataPtr(cctkGH, 0, "CT_Analytic::testdzK");
plist[3] = (CCTK_REAL*) CCTK_VarDataPtr(cctkGH, 0, "CT_MultiLevel::ct_a0[0]");
plist[4] = (CCTK_REAL*) CCTK_VarDataPtr(cctkGH, 0, "CT_Analytic::testc0");
