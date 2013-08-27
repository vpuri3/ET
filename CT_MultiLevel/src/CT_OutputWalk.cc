#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "carpet.hh"
#include "loopcontrol.h"
#include "Carpet/Carpet/src/modes.hh"

#include "CT_MultiLevel.hh"

extern "C" void CT_OutputWalk(CCTK_ARGUMENTS)
{
  BEGIN_MAP_LOOP(cctkGH,CCTK_GF) {
    BEGIN_COMPONENT_LOOP(cctkGH, CCTK_GF) {
      DECLARE_CCTK_ARGUMENTS;
      DECLARE_CCTK_PARAMETERS;

      if (CT_ProcessOwnsData())
      {
        CCTK_REAL c_value, r_value;
        CCTK_INT npoints = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];

        for (int nequation=0; nequation < number_of_equations; nequation++)
        {
          LC_LOOP3 (EL_INI,
                    i, j, k, 0, 0, 0, cctk_lsh[0], cctk_lsh[1], cctk_lsh[2], 
                    cctk_lsh[0], cctk_lsh[1], cctk_lsh[2])
          {
            int sindex = CCTK_GFINDEX3D(cctkGH,i,j,k);
            int index = sindex + nequation * npoints;

	    if (r[sindex] == 0)
            {
              c_value = ct_c1[index];
              r_value = ct_residual[index];
            }
          } LC_ENDLOOP3 (EL_INI);

          char fname[100];
          sprintf(fname, "walk_eqn%d.asc", nequation);
          CT_WritePairs(c_value, r_value, fname); 
        }
      }

    } END_COMPONENT_LOOP;
  } END_MAP_LOOP;

  return;
}
