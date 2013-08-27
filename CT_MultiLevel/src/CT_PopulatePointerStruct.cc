#include <string.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "carpet.hh"
#include "loopcontrol.h"
#include "Carpet/Carpet/src/modes.hh"

#include "CT_MultiLevel.hh"

extern "C" void CT_PopulatePointerStruct(CCTK_ARGUMENTS, struct coeffptr *cptr)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  CCTK_REAL **plist = (CCTK_REAL**) malloc(sizeof(CCTK_REAL*));

  if(CCTK_Equals(model, "Bowen-York"))
  {
    #include "extras/BowenYork.cc"
  }
  else if (CCTK_Equals(model, "Expanding BH lattice"))
  {
    #include "extras/ExpLat.cc"
  }
  else if (CCTK_Equals(model, "Lump"))
  {
    #include "extras/Lump.cc"
  }

  for (int nequation=0; nequation < number_of_equations; nequation++)
  {
    struct coeffptr tmp = {
      (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, cxx_gfname[nequation]),
      (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, cxy_gfname[nequation]),
      (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, cxz_gfname[nequation]),
      (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, cyy_gfname[nequation]),
      (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, cyz_gfname[nequation]),
      (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, czz_gfname[nequation]),
      (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, cx_gfname[nequation]),
      (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, cy_gfname[nequation]),
      (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, cz_gfname[nequation]),
      (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, c0_gfname[nequation]),
      (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, c1_gfname[nequation]),
      (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, c2_gfname[nequation]),
      (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, c3_gfname[nequation]),
      (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, c4_gfname[nequation]),
      (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, a0_gfname[nequation]),
      (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, a1_gfname[nequation]),
      (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, a2_gfname[nequation]),
      (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, a3_gfname[nequation]),
      (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, a4_gfname[nequation]),
      ct_cxx,
      ct_cxy,
      ct_cxz,
      ct_cyy,
      ct_cyz,
      ct_czz,
      ct_cx,
      ct_cy,
      ct_cz,
      ct_c0,
      ct_c1,
      ct_c2,
      ct_c3,
      ct_c4,
      ct_a0,
      ct_a1,
      ct_a2,
      ct_a3,
      ct_a4,
      ct_psi,
      ct_auxiliary,
      plist 
    };
    cptr[nequation] = tmp;
  }

  return;
}
