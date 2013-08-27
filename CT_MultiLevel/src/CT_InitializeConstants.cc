#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "carpet.hh"
#include "loopcontrol.h"
#include "Carpet/Carpet/src/modes.hh"

#include "CT_MultiLevel.hh"

extern "C" void CT_InitializeConstants(CCTK_ARGUMENTS)
{
  BEGIN_REFLEVEL_LOOP(cctkGH) {
    BEGIN_MAP_LOOP(cctkGH,CCTK_GF) {
      BEGIN_COMPONENT_LOOP(cctkGH, CCTK_GF) {
        DECLARE_CCTK_ARGUMENTS;
        DECLARE_CCTK_PARAMETERS;

        if (CT_ProcessOwnsData())
        {
          LC_LOOP3 (EL_INI,
                    i, j, k, 0, 0, 0, cctk_lsh[0], cctk_lsh[1], cctk_lsh[2], 
                    cctk_lsh[0], cctk_lsh[1], cctk_lsh[2])
          {
            int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

            ct_zero[index] = 0;
          } LC_ENDLOOP3 (EL_INI);
        }

      } END_COMPONENT_LOOP;
    } END_MAP_LOOP;
  } END_REFLEVEL_LOOP; 

  return;
}
