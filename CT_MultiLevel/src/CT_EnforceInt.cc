#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "carpet.hh"
#include "loopcontrol.h"
#include "Carpet/Carpet/src/modes.hh"

#include "CT_MultiLevel.hh"

extern "C" void CT_EnforceInt(CCTK_ARGUMENTS, CCTK_REAL *value)
{
  BEGIN_MAP_LOOP(cctkGH,CCTK_GF) {
    BEGIN_COMPONENT_LOOP(cctkGH, CCTK_GF) {
      DECLARE_CCTK_ARGUMENTS;
      DECLARE_CCTK_PARAMETERS;

      if (CT_ProcessOwnsData())
      {
        CCTK_INT npoints = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];

        if(CCTK_Equals(model, "Lump"))
        {
          #include "integral/Lump.cc"
        } else if(CCTK_Equals(model, "Expanding BH lattice")) {
          if (CCTK_Equals(reset_psi, "to value"))
          {
            #include "integral/ExpLat.cc"
          }
          if (CCTK_Equals(reset_psi, "through integrability"))
          {
            #include "integral/ExpLatNewt.cc"
          }
        } else if(CCTK_Equals(model, "Inhomogeneous Helmholtz")) {
          if (CCTK_Equals(reset_psi, "to value"))
          {
            #include "integral/IHelmholtzCoeff.cc"
          }
          if (CCTK_Equals(reset_psi, "through integrability"))
          {
            #include "integral/IHelmholtz.cc"
          }
        }
      }

    } END_COMPONENT_LOOP;
  } END_MAP_LOOP;

  return;
}
