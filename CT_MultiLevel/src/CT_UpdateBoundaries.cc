#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include "cctk_Schedule.h"

#include "carpet.hh"
#include "loopcontrol.h"
#include "Carpet/Carpet/src/modes.hh"

#include "CT_MultiLevel.hh"

extern "C" void CT_UpdateBoundaries(CCTK_ARGUMENTS, const char *varname)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_Equals(veryverbose, "yes")) CCTK_VInfo(CCTK_THORNSTRING, "Updating boundaries of %s...", varname);
   
//  BEGIN_MAP_LOOP(cctkGH,CCTK_GF) {
//    BEGIN_COMPONENT_LOOP(cctkGH, CCTK_GF) {
//      if (CT_ProcessOwnsData())
//      {
        if (CCTK_SyncGroup(cctkGH, varname) < 0)
          CCTK_VWarn(0, __LINE__, __FILE__, "Couldn't sync variable %s", varname); 
//      }
//    } END_COMPONENT_LOOP;
//  } END_MAP_LOOP;

  if (!CCTK_Equals(boundary_conditions, "none") && CCTK_Equals(varname, "CT_MultiLevel::psi"))
  {
    CT_Boundaries(CCTK_PASS_CTOC, varname);
    //CCTK_SyncGroup(cctkGH, varname) // Need this?    
  } else {
    int ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, cctk_nghostzones[0], -1, varname, "none");
    if (ierr < 0)
      CCTK_VWarn(0, __LINE__, __FILE__, "Failed to register boundary conditions for %s.", varname);

    CCTK_ScheduleTraverse("ApplyBCs", cctkGH, NULL);
  }

  return;
}
