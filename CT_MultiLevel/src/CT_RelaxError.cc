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

extern "C" void CT_RelaxError(CCTK_ARGUMENTS, 
                              CCTK_INT level, 
                              CCTK_INT nsteps,
                              CCTK_INT init_err, 
                              CCTK_INT prol_psi, 
                              CCTK_INT rest_psi,
                              CCTK_INT rest_res,
                              CCTK_INT add_err,
                              CCTK_INT rset_psi,
                              CCTK_INT init_coeff)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_VInfo(CCTK_THORNSTRING, "=== Level %d ===", level);

  ENTER_LEVEL_MODE(cctkGH,level) {
    CCTK_INT step = 0;
    CCTK_REAL norm = 1e+03, tmpnorm;
    CCTK_REAL *zero = (CCTK_REAL *) malloc(number_of_equations*sizeof(CCTK_REAL));

//    if (init_coeff) CT_InitializeCoefficients(CCTK_PASS_CTOC);

    if (prol_psi) CT_Prolongate(CCTK_PASS_CTOC, "CT_MultiLevel::psi");
    if (rest_psi)
    {
      CT_Restrict(CCTK_PASS_CTOC, "CT_MultiLevel::psi"); 
      CT_Restrict(CCTK_PASS_CTOC, "CT_MultiLevel::coeffs"); 
      CT_Restrict(CCTK_PASS_CTOC, "CT_Analytic::CT_testK"); 
      CT_Restrict(CCTK_PASS_CTOC, "CT_Analytic::CT_testdxK"); 
      CT_Restrict(CCTK_PASS_CTOC, "CT_Analytic::CT_testdyK"); 
      CT_Restrict(CCTK_PASS_CTOC, "CT_Analytic::CT_testdzK"); 
      CT_Restrict(CCTK_PASS_CTOC, "CT_Analytic::CT_testc1"); 
      CT_Restrict(CCTK_PASS_CTOC, "CT_Analytic::CT_testc2"); 
    }
    if (rest_res)
    {
      CT_Restrict(CCTK_PASS_CTOC, "CT_MultiLevel::residual"); 
      for (int nequation=0; nequation < number_of_equations; nequation++)
        CT_CopyResidual(CCTK_PASS_CTOC, nequation);
    }
    if (add_err < 0)
    {
      CT_Restrict(CCTK_PASS_CTOC, "CT_MultiLevel::err"); 
      CT_AddErrorToPsi(CCTK_PASS_CTOC);
    }
    else if (add_err > 0)
    {
      CT_Prolongate(CCTK_PASS_CTOC, "CT_MultiLevel::err"); 
      CT_AddErrorToPsi(CCTK_PASS_CTOC);
    }
    if (init_err) CT_InitializeError(CCTK_PASS_CTOC);
    else CT_RestoreError(CCTK_PASS_CTOC);
    CT_UpdateBoundaries(CCTK_PASS_CTOC, "CT_MultiLevel::psi");
    for (int nequation=0; nequation < number_of_equations; nequation++)
      CT_Copy(CCTK_PASS_CTOC, "CT_MultiLevel::ct_residual_above[0]", "CT_MultiLevel::ct_residual_above_copy[0]", nequation);

    while (norm > tol && step < nsteps)
    {
      //if (CCTK_Equals(disable[nequation], "no"))
      CT_SolveErrorEquation(CCTK_PASS_CTOC, &norm, &step);
      for (int nequation=0; nequation < number_of_equations; nequation++) zero[nequation] = 0;

      if (rset_psi) CT_Reset(CCTK_PASS_CTOC, "CT_MultiLevel::ct_err", zero, 0);
      CT_UpdateBoundaries(CCTK_PASS_CTOC, "CT_MultiLevel::err");

//      CT_InitializeCoefficients(CCTK_PASS_CTOC);

      if (CCTK_Equals(compare_to_exact, "yes")) CT_CompareToExact(CCTK_PASS_CTOC, 1);
      if (CCTK_Equals(output_walk,"yes"))
      {
        //CT_CalcPsiResidual(CCTK_PASS_CTOC, step, 0);
        CT_OutputWalk(CCTK_PASS_CTOC);
      }
      CT_CalcErrResidual(CCTK_PASS_CTOC, step, 0, &norm);
    }

    for (int nequation=0; nequation < number_of_equations; nequation++)
      CT_Copy(CCTK_PASS_CTOC, "CT_MultiLevel::ct_err[0]", "CT_MultiLevel::ct_err_copy[0]", nequation);
    
    double tmp;
    CT_CalcErrResidual(CCTK_PASS_CTOC, step, 1, &tmp);
    CT_UpdateBoundaries(CCTK_PASS_CTOC, "CT_MultiLevel::residual");
    CT_UpdateBoundaries(CCTK_PASS_CTOC, "CT_MultiLevel::coeffs");
    for (int nequation=0; nequation < number_of_equations; nequation++)
    {
      CT_Copy(CCTK_PASS_CTOC, "CT_MultiLevel::ct_residual[0]", "CT_MultiLevel::ct_residual_copy[0]", nequation);

      /*if (add_err_psi)
      {
        CT_AddErrorToPsi(CCTK_PASS_CTOC);
        if (rset_psi) CT_Reset(CCTK_PASS_CTOC, CCTK_VarIndex("CT_MultiLevel::ct_psi"), 1, 1);
      }*/
      CT_Copy(CCTK_PASS_CTOC, "CT_MultiLevel::ct_psi[0]", "CT_MultiLevel::ct_psi_copy[0]", nequation);
    }

  } LEAVE_LEVEL_MODE; 

  return;
}
