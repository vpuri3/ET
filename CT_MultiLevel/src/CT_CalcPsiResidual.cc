#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "carpet.hh"
#include "loopcontrol.h"
#include "Carpet/Carpet/src/modes.hh"

#include "CT_MultiLevel.hh"

extern "C" void CT_CalcPsiResidual(CCTK_ARGUMENTS, CCTK_INT step, CCTK_INT output, CCTK_REAL *norm)
{
  DECLARE_CCTK_PARAMETERS;

  CCTK_REAL eqnnorm;
  *norm = 0;

  for (int nequation=0; nequation < number_of_equations; nequation++)
  {
    BEGIN_MAP_LOOP(cctkGH,CCTK_GF) {
      BEGIN_COMPONENT_LOOP(cctkGH, CCTK_GF) {
        DECLARE_CCTK_ARGUMENTS;

        if (CCTK_Equals(verbose,"yes") && output) CCTK_Info(CCTK_THORNSTRING, "Calculating the residual for psi");

        if (CT_ProcessOwnsData())
        {
          int xgh = cctk_nghostzones[0];
          int ygh = cctk_nghostzones[1];
          int zgh = cctk_nghostzones[2];
          int imin = xgh; 
          int imax = cctk_lsh[0]-xgh;
          int jmin = ygh;
          int jmax = cctk_lsh[1]-ygh;
          int kmin = zgh;
          int kmax = cctk_lsh[2]-zgh;

          CCTK_INT npoints = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];

          struct coeffptr *cptr = (struct coeffptr *) malloc(number_of_equations*sizeof(struct coeffptr));
          CT_PopulatePointerStruct(CCTK_PASS_CTOC, cptr);

          LC_LOOP3 (EL_PRS,
                    i, j, k, imin, jmin, kmin, imax, jmax, kmax,
                    cctk_lsh[0], cctk_lsh[1], cctk_lsh[2])
          {
            int sindex = CCTK_GFINDEX3D(cctkGH,i,j,k);
            int index = CCTK_GFINDEX3D(cctkGH,i,j,k) + nequation*npoints;
            struct stencil stnc = {
              (int)CCTK_GFINDEX3D(cctkGH,i,j,k),
              (int)CCTK_GFINDEX3D(cctkGH,i-1,j,k),
              (int)CCTK_GFINDEX3D(cctkGH,i+1,j,k),
              (int)CCTK_GFINDEX3D(cctkGH,i,j-1,k),
              (int)CCTK_GFINDEX3D(cctkGH,i,j+1,k),
              (int)CCTK_GFINDEX3D(cctkGH,i,j,k-1),
              (int)CCTK_GFINDEX3D(cctkGH,i,j,k+1),
              (int)CCTK_GFINDEX3D(cctkGH,i-2,j,k),
              (int)CCTK_GFINDEX3D(cctkGH,i+2,j,k),
              (int)CCTK_GFINDEX3D(cctkGH,i,j-2,k),
              (int)CCTK_GFINDEX3D(cctkGH,i,j+2,k),
              (int)CCTK_GFINDEX3D(cctkGH,i,j,k-2),
              (int)CCTK_GFINDEX3D(cctkGH,i,j,k+2),
              (int)CCTK_GFINDEX3D(cctkGH,i-1,j-1,k),
              (int)CCTK_GFINDEX3D(cctkGH,i-1,j+1,k),
              (int)CCTK_GFINDEX3D(cctkGH,i+1,j-1,k),
              (int)CCTK_GFINDEX3D(cctkGH,i+1,j+1,k),
              (int)CCTK_GFINDEX3D(cctkGH,i-1,j,k-1),
              (int)CCTK_GFINDEX3D(cctkGH,i-1,j,k+1),
              (int)CCTK_GFINDEX3D(cctkGH,i+1,j,k-1),
              (int)CCTK_GFINDEX3D(cctkGH,i+1,j,k+1),
              (int)CCTK_GFINDEX3D(cctkGH,i,j-1,k-1),
              (int)CCTK_GFINDEX3D(cctkGH,i,j-1,k+1),
              (int)CCTK_GFINDEX3D(cctkGH,i,j+1,k-1),
              (int)CCTK_GFINDEX3D(cctkGH,i,j+1,k+1),
              CCTK_DELTA_SPACE(0),
              CCTK_DELTA_SPACE(1),
              CCTK_DELTA_SPACE(2),
              fd_order,
              npoints
            };
            CT_InitializeCoefficients(sindex, nequation, &cptr[nequation], stnc);

            CCTK_REAL psix, psiy, psiz, psixx, psiyy, psizz, psixy, psixz, psiyz;
            CT_FD(ct_psi, stnc, nequation,
                  &psix, &psiy, &psiz, &psixx, &psiyy, &psizz, &psixy, &psixz, &psiyz);

            CCTK_REAL psi0 = pow(ct_psi[index]+ct_a0[index], n0[nequation]);
            CCTK_REAL psi1 = pow(ct_psi[index]+ct_a1[index], n1[nequation]);
            CCTK_REAL psi2 = pow(ct_psi[index]+ct_a2[index], n2[nequation]);
            CCTK_REAL psi3 = pow(ct_psi[index]+ct_a3[index], n3[nequation]);
            CCTK_REAL psi4 = pow(ct_psi[index]+ct_a4[index], n4[nequation]);

            ct_residual[index] = - ct_cxx[index] * psixx
                                 - ct_cyy[index] * psiyy
                                 - ct_czz[index] * psizz
                                 - ct_cxy[index] * psixy
                                 - ct_cxz[index] * psixz
                                 - ct_cyz[index] * psiyz
                                 - ct_cx[index] * psix
                                 - ct_cy[index] * psiy 
                                 - ct_cz[index] * psiz
                                 - ct_c0[index] * psi0
                                 - ct_c1[index] * psi1
                                 - ct_c2[index] * psi2
                                 - ct_c3[index] * psi3
                                 - ct_c4[index] * psi4;
          } LC_ENDLOOP3 (EL_PRS);
          free(cptr);
        } // if
      } END_COMPONENT_LOOP;
    } END_MAP_LOOP;

    CT_Norm(CCTK_PASS_CTOC, "CT_MultiLevel::ct_residual[0]", &eqnnorm, nequation);

    if (output) CCTK_VInfo(CCTK_THORNSTRING, " * Equation #%d: %d iterations, norm of final residual = %1.10e", nequation, step, eqnnorm);
    if (eqnnorm > *norm) *norm = eqnnorm;
  } // for
  return;
}
