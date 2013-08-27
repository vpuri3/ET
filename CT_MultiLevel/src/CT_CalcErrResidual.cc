#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "carpet.hh"
#include "loopcontrol.h"
#include "Carpet/Carpet/src/modes.hh"

#include "CT_MultiLevel.hh"

extern "C" void CT_CalcErrResidual(CCTK_ARGUMENTS, CCTK_INT step)
{
  DECLARE_CCTK_PARAMETERS;

  for (int nequation=0; nequation < number_of_equations; nequation++)
  {
    BEGIN_MAP_LOOP(cctkGH,CCTK_GF) {
      BEGIN_COMPONENT_LOOP(cctkGH, CCTK_GF) {
        DECLARE_CCTK_ARGUMENTS;

        if (CCTK_Equals(verbose,"yes")) CCTK_Info(CCTK_THORNSTRING, "Calculating the residual for err");

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

          LC_LOOP3 (EL_ERS,
                    i, j, k, imin, jmin, kmin, imax, jmax, kmax,
                    cctk_lsh[0], cctk_lsh[1], cctk_lsh[2])
          {
            int sindex = CCTK_GFINDEX3D(cctkGH,i,j,k);
            int index = CCTK_GFINDEX3D(cctkGH,i,j,k) + nequation*npoints;
            struct stencil stnc = {
              CCTK_GFINDEX3D(cctkGH,i,j,k),
              CCTK_GFINDEX3D(cctkGH,i-1,j,k),
              CCTK_GFINDEX3D(cctkGH,i+1,j,k),
              CCTK_GFINDEX3D(cctkGH,i,j-1,k),
              CCTK_GFINDEX3D(cctkGH,i,j+1,k),
              CCTK_GFINDEX3D(cctkGH,i,j,k-1),
              CCTK_GFINDEX3D(cctkGH,i,j,k+1),
              CCTK_GFINDEX3D(cctkGH,i-2,j,k),
              CCTK_GFINDEX3D(cctkGH,i+2,j,k),
              CCTK_GFINDEX3D(cctkGH,i,j-2,k),
              CCTK_GFINDEX3D(cctkGH,i,j+2,k),
              CCTK_GFINDEX3D(cctkGH,i,j,k-2),
              CCTK_GFINDEX3D(cctkGH,i,j,k+2),
              CCTK_GFINDEX3D(cctkGH,i-1,j-1,k),
              CCTK_GFINDEX3D(cctkGH,i-1,j+1,k),
              CCTK_GFINDEX3D(cctkGH,i+1,j-1,k),
              CCTK_GFINDEX3D(cctkGH,i+1,j+1,k),
              CCTK_GFINDEX3D(cctkGH,i-1,j,k-1),
              CCTK_GFINDEX3D(cctkGH,i-1,j,k+1),
              CCTK_GFINDEX3D(cctkGH,i+1,j,k-1),
              CCTK_GFINDEX3D(cctkGH,i+1,j,k+1),
              CCTK_GFINDEX3D(cctkGH,i,j-1,k-1),
              CCTK_GFINDEX3D(cctkGH,i,j-1,k+1),
              CCTK_GFINDEX3D(cctkGH,i,j+1,k-1),
              CCTK_GFINDEX3D(cctkGH,i,j+1,k+1),
              CCTK_DELTA_SPACE(0),
              CCTK_DELTA_SPACE(1),
              CCTK_DELTA_SPACE(2),
              fd_order,
              npoints
            };
            CT_InitializeCoefficients(sindex, nequation, &cptr[nequation], stnc);

            CCTK_REAL errx, erry, errz, errxx, erryy, errzz, errxy, errxz, erryz;
            CT_FD(ct_err, stnc, nequation,
                  &errx, &erry, &errz, &errxx, &erryy, &errzz, &errxy, &errxz, &erryz);

            CCTK_REAL err0 = pow(ct_psi[index] + ct_a0[index] + ct_err[index], n0[nequation]) - pow(ct_psi[index] + ct_a0[index], n0[nequation]);
            CCTK_REAL err1 = pow(ct_psi[index] + ct_a1[index] + ct_err[index], n1[nequation]) - pow(ct_psi[index] + ct_a1[index], n1[nequation]);
            CCTK_REAL err2 = pow(ct_psi[index] + ct_a2[index] + ct_err[index], n2[nequation]) - pow(ct_psi[index] + ct_a2[index], n2[nequation]);
            CCTK_REAL err3 = pow(ct_psi[index] + ct_a3[index] + ct_err[index], n3[nequation]) - pow(ct_psi[index] + ct_a3[index], n3[nequation]);
            CCTK_REAL err4 = pow(ct_psi[index] + ct_a4[index] + ct_err[index], n4[nequation]) - pow(ct_psi[index] + ct_a4[index], n4[nequation]);
            CCTK_REAL res  = ct_residual_above[index];

            ct_residual[index] = - ct_cxx[index] * errxx
                                 - ct_cyy[index] * erryy
                                 - ct_czz[index] * errzz
                                 - ct_cxy[index] * errxy
                                 - ct_cxz[index] * errxz
                                 - ct_cyz[index] * erryz
                                 - ct_cx[index] * errx
                                 - ct_cy[index] * erry 
                                 - ct_cz[index] * errz
                                 - ct_c0[index] * err0
                                 - ct_c1[index] * err1
                                 - ct_c2[index] * err2
                                 - ct_c3[index] * err3
                                 - ct_c4[index] * err4
                                 + res;
          } LC_ENDLOOP3 (EL_ERS);
          free(cptr);
        } // if
      } END_COMPONENT_LOOP;
    } END_MAP_LOOP;

    CCTK_REAL norm;
    CT_Norm(CCTK_PASS_CTOC, "CT_MultiLevel::ct_residual[0]", &norm, nequation);

    CCTK_VInfo(CCTK_THORNSTRING, " * Equation #%d: %d iterations, norm of final residual = %1.10e", nequation, step, norm);
  } // for

  return;
}
