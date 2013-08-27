#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "carpet.hh"
#include "loopcontrol.h"
#include "Carpet/Carpet/src/modes.hh"

#include "CT_MultiLevel.hh"

extern "C" void CT_RestoreError(CCTK_ARGUMENTS)
{
  BEGIN_MAP_LOOP(cctkGH,CCTK_GF) {
    BEGIN_COMPONENT_LOOP(cctkGH, CCTK_GF) {
      DECLARE_CCTK_ARGUMENTS;
      DECLARE_CCTK_PARAMETERS;

      if (CCTK_Equals(verbose, "yes")) CCTK_Info(CCTK_THORNSTRING, "Restoring error...");

      if (CT_ProcessOwnsData())
      {
        for (int nequation=0; nequation < number_of_equations; nequation++)
        {
          CCTK_INT npoints = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];

          LC_LOOP3 (EL_INI,
                    i, j, k, 0, 0, 0, cctk_lsh[0], cctk_lsh[1], cctk_lsh[2], 
                    cctk_lsh[0], cctk_lsh[1], cctk_lsh[2])
          {
            int index = CCTK_GFINDEX3D(cctkGH,i,j,k) + nequation * npoints;

            ct_err[index] = ct_err_copy[index];
          } LC_ENDLOOP3 (EL_INI);
        }
      }

    } END_COMPONENT_LOOP;
  } END_MAP_LOOP;

  return;
}
