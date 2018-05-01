#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include "cctk_Schedule.h"

#include "carpet.hh"
#include "loopcontrol.h"
#include "Carpet/Carpet/src/modes.hh"

#include "CT_MultiLevel.hh"

using namespace Carpet;

extern "C" void CT_V(CCTK_ARGUMENTS, CCTK_INT toplevel, CCTK_INT init_psi)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT nlevels = 2*(toplevel+1)-1;
  CCTK_INT rset_psi = CCTK_Equals(reset_psi, "no")? 0 : CCTK_Equals(reset_psi, "to value")? 1 : -1;

  CCTK_INT *levels   = (CCTK_INT*)malloc(nlevels*sizeof(CCTK_INT));
  CCTK_INT *downward = (CCTK_INT*)malloc(nlevels*sizeof(CCTK_INT));
  CT_FillFlagArrays(levels, downward, toplevel);

  for (int irl=0; irl<nlevels; irl++)
  {
    CCTK_INT rl = levels[irl];

    CCTK_INT nsteps = rl == 0 ? nrelsteps_bottom : rl == toplevel ? nrelsteps_top : downward[irl] ? nrelsteps_down : nrelsteps_up;

    if (downward[irl])
    {
      if (rl == toplevel)
        // FIXME: here we're using InitializePsi even when the 
        // prolongation of psi occurs, just to set the boundary
        // conditions before the prolongation happens. Ideally
        // the initialization of interior and boundary points
        // should happen separately.
        CT_RelaxPsi(CCTK_PASS_CTOC, 
                    rl, 
                    nsteps, 
                    1,            // Initialize psi
                    !init_psi,    // Prolongate psi
                    0,            // Restrict psi
                    0,            // Add the error to psi
                    rset_psi,     // Reset psi
                    1,            // Initialize the coefficients
                    enforce_int); // Enforce the integral constraint
      else
        CT_RelaxError(CCTK_PASS_CTOC, 
                      rl, 
                      nsteps, 
                      1,          // Initialize the error 
                      0,          // Prolongate psi
                      1,          // Restrict psi
                      1,          // Restrict the residual
                      -(rl!=toplevel-1),  // Add the error to psi
                      rset_psi,   // Reset the error
                      !enforce_int); // Initialize the coefficients
    } 
    else
    {
      if (rl != toplevel)
        //CT_RelaxError(CCTK_PASS_CTOC, rl, nsteps, 0, 0, 0, 0, (rl>0), rset_psi, 0);
        CT_RelaxError(CCTK_PASS_CTOC, 
                      rl, 
                      nsteps, 
                      0,            // Initialize the error
                      0,            // Prolongate psi
                      0,            // Restrict psi
                      0,            // Restrict the residual
                      (rl>0),       // Add the error to psi 
                      rset_psi,     // Reset psi 
                      1);           // Initialize the coefficients
      else
        CT_RelaxPsi(CCTK_PASS_CTOC, 
                    rl, 
                    nsteps, 
                    0,              // Initialize psi 
                    0,              // Prolongate psi
                    0,              // Restrict psi
                    1,              // Add error to psi
                    rset_psi,       // Reset psi
                    1,              // Initialize the coefficients
                    enforce_int);   // Enforce the integral constraint
    }
  } 

  free(levels);
  free(downward);

  return;
}
