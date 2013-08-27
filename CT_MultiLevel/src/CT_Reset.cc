#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "carpet.hh"
#include "loopcontrol.h"
#include "Carpet/Carpet/src/modes.hh"

#include "CT_MultiLevel.hh"

extern "C" void CT_Reset(CCTK_ARGUMENTS, char *varName, CCTK_REAL *resetvalue, CCTK_INT relative)
{
  BEGIN_MAP_LOOP(cctkGH,CCTK_GF) {
    BEGIN_COMPONENT_LOOP(cctkGH, CCTK_GF) {
      DECLARE_CCTK_ARGUMENTS;
      DECLARE_CCTK_PARAMETERS;

      if (CT_ProcessOwnsData())
      {
        for (int nequation=0; nequation < number_of_equations; nequation++)
        {
          CCTK_REAL value;
          CCTK_INT npoints = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];

          char varName1[100];
          sprintf(varName1, "%s[0]", varName);
          CCTK_REAL *gf = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, varName1);
          sprintf(varName1, "%s[%d]", varName, nequation);
          if (relative)
          {
            value = 0;
            if (CCTK_Equals(veryverbose, "yes"))
              CCTK_VInfo(CCTK_THORNSTRING, "Resetting the value of %s by %f...", varName, resetvalue[nequation]);
          }
          else
          {
            value = CT_GetValue(CCTK_PASS_CTOC, reset_x, reset_y, reset_z, varName1);
            if (CCTK_Equals(veryverbose, "yes"))
              CCTK_VInfo(CCTK_THORNSTRING, "Resetting the value of %s at (%f, %f, %f) to %f...", varName, reset_x, reset_y, reset_z, resetvalue[nequation]);
          }

          LC_LOOP3 (EL_INI,
                    i, j, k, 0, 0, 0, cctk_lsh[0], cctk_lsh[1], cctk_lsh[2], 
                    cctk_lsh[0], cctk_lsh[1], cctk_lsh[2])
          {
            int index = CCTK_GFINDEX3D(cctkGH,i,j,k) + nequation * npoints;

            gf[index] = gf[index] - value + resetvalue[nequation];
          } LC_ENDLOOP3 (EL_INI);
        }
      }

    } END_COMPONENT_LOOP;
  } END_MAP_LOOP;

  return;
}
