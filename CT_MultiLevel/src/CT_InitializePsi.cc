#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "carpet.hh"
#include "loopcontrol.h"
#include "Carpet/Carpet/src/modes.hh"

#include "CT_MultiLevel.hh"

extern "C" void CT_InitializePsi(CCTK_ARGUMENTS)
{
  BEGIN_MAP_LOOP(cctkGH,CCTK_GF) {
    BEGIN_COMPONENT_LOOP(cctkGH, CCTK_GF) {
      DECLARE_CCTK_ARGUMENTS;
      DECLARE_CCTK_PARAMETERS;

      if (CT_ProcessOwnsData())
      {
        if (CCTK_Equals(verbose, "yes")) CCTK_Info(CCTK_THORNSTRING, "Initializing psi...");

        CCTK_INT npoints = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];

        for (int nequation=0; nequation < number_of_equations; nequation++)
        {
          CCTK_REAL *pvar = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, inipsi_gfname[nequation]);

          LC_LOOP3 (EL_INI,
                    i, j, k, 0, 0, 0, cctk_lsh[0], cctk_lsh[1], cctk_lsh[2], 
                    cctk_lsh[0], cctk_lsh[1], cctk_lsh[2])
          {
            int sindex = CCTK_GFINDEX3D(cctkGH,i,j,k);
            int index = sindex + nequation * npoints;

            ct_psi[index] = pvar[sindex];

          } LC_ENDLOOP3 (EL_INI);
        }
      }

    } END_COMPONENT_LOOP;
  } END_MAP_LOOP;

  return;
}
