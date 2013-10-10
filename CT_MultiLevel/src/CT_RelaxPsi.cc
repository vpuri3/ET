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

extern "C" void CT_RelaxPsi(CCTK_ARGUMENTS, 
                            CCTK_INT level, 
                            CCTK_INT nsteps,
                            CCTK_INT init_psi, 
                            CCTK_INT prol_psi, 
                            CCTK_INT rest_psi,
                            CCTK_INT add_err,
                            CCTK_INT rset_psi,
                            CCTK_INT init_coeff,
                            CCTK_INT enf_int)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_VInfo(CCTK_THORNSTRING, "=== Level %d ===", level);

  ENTER_LEVEL_MODE(cctkGH,level) {
    CCTK_INT step = 0;
    CCTK_REAL norm = 1e+03, tmpnorm;
    CCTK_REAL *value = (CCTK_REAL *) malloc(number_of_equations*sizeof(CCTK_REAL));

    if (init_psi) CT_InitializePsi(CCTK_PASS_CTOC);
    if (prol_psi) CT_Prolongate(CCTK_PASS_CTOC, "CT_MultiLevel::psi"); 
    if (rest_psi)
    {
      CT_Restrict(CCTK_PASS_CTOC, "CT_MultiLevel::psi"); 
      CT_Restrict(CCTK_PASS_CTOC, "CT_MultiLevel::coeffs"); 
      CT_Restrict(CCTK_PASS_CTOC, "CT_Analytic::CT_testK"); 
      CT_Restrict(CCTK_PASS_CTOC, "CT_Analytic::CT_testdxK"); 
      CT_Restrict(CCTK_PASS_CTOC, "CT_Analytic::CT_testdyK"); 
      CT_Restrict(CCTK_PASS_CTOC, "CT_Analytic::CT_testdzK"); 
    }
    if (add_err)
    {
      CT_InitializeError(CCTK_PASS_CTOC);
      CT_Prolongate(CCTK_PASS_CTOC, "CT_MultiLevel::err"); 
      CT_AddErrorToPsi(CCTK_PASS_CTOC);
    }
//    if (init_coeff) CT_InitializeCoefficients(CCTK_PASS_CTOC);

    // In this case, unlike the CT_RelaxError case, the norm shouldn't be 
    // used as a stopping criterion, but rather the norm of the residual.
    while (norm > tol && step < nsteps)
    {
      CT_SolvePsiEquation(CCTK_PASS_CTOC, &norm, &step);

      if (rset_psi > 0 && step%reset_every == 0)
      {
        for (int nequation=0; nequation < number_of_equations; nequation++) value[nequation] = reset_value[nequation];
        CT_Reset(CCTK_PASS_CTOC, "CT_MultiLevel::ct_psi", value, 0);
      }
      if (enf_int) CT_EnforceInt(CCTK_PASS_CTOC, value); // Is the conditional necessary?
      if (rset_psi < 0 && step%reset_every == 0)
      {
        CT_Reset(CCTK_PASS_CTOC, "CT_MultiLevel::ct_psi", value, 1);
      }
      CT_UpdateBoundaries(CCTK_PASS_CTOC, "CT_MultiLevel::psi");

      if (CCTK_Equals(compare_to_exact, "yes")) CT_CompareToExact(CCTK_PASS_CTOC, 0);
      if (CCTK_Equals(output_walk,"yes"))
      {
        CT_CalcPsiResidual(CCTK_PASS_CTOC, step, 0);
        CT_OutputWalk(CCTK_PASS_CTOC);
      }
    }

    for (int nequation=0; nequation < number_of_equations; nequation++)
      CT_Copy(CCTK_PASS_CTOC, "CT_MultiLevel::ct_psi[0]", "CT_MultiLevel::ct_psi_copy[0]", nequation);

    CT_CalcPsiResidual(CCTK_PASS_CTOC, step, 1);
    CT_UpdateBoundaries(CCTK_PASS_CTOC, "CT_MultiLevel::residual");
    CT_UpdateBoundaries(CCTK_PASS_CTOC, "CT_MultiLevel::coeffs");
    for (int nequation=0; nequation < number_of_equations; nequation++)
      CT_Copy(CCTK_PASS_CTOC, "CT_MultiLevel::ct_residual[0]", "CT_MultiLevel::ct_residual_copy[0]", nequation);

  } LEAVE_LEVEL_MODE; 

  return;
}
