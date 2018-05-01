plist = (CCTK_REAL**) realloc(plist, 6*sizeof(CCTK_REAL*));

plist[0] = (CCTK_REAL*) CCTK_VarDataPtr(cctkGH, 0, "CT_Analytic::testAxx");
plist[1] = (CCTK_REAL*) CCTK_VarDataPtr(cctkGH, 0, "CT_Analytic::testAxy");
plist[2] = (CCTK_REAL*) CCTK_VarDataPtr(cctkGH, 0, "CT_Analytic::testAxz");
plist[3] = (CCTK_REAL*) CCTK_VarDataPtr(cctkGH, 0, "CT_Analytic::testAyy");
plist[4] = (CCTK_REAL*) CCTK_VarDataPtr(cctkGH, 0, "CT_Analytic::testAyz");
plist[5] = (CCTK_REAL*) CCTK_VarDataPtr(cctkGH, 0, "CT_Analytic::testAzz");
