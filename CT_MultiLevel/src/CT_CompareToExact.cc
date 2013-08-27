#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "carpet.hh"
#include "loopcontrol.h"
#include "Carpet/Carpet/src/modes.hh"

#include "CT_MultiLevel.hh"

extern "C" void CT_CompareToExact(CCTK_ARGUMENTS, CCTK_INT outdated_psi)
{
  DECLARE_CCTK_PARAMETERS;

  for (int nequation=0; nequation < number_of_equations; nequation++)
  {
    BEGIN_MAP_LOOP(cctkGH,CCTK_GF) {
      BEGIN_COMPONENT_LOOP(cctkGH, CCTK_GF) {
        DECLARE_CCTK_ARGUMENTS;

        if (CT_ProcessOwnsData())
        {
          if (CCTK_Equals(veryverbose, "yes")) CCTK_Info(CCTK_THORNSTRING, "Comparing with the exact solution...");

          CCTK_REAL *epsi = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, exact_solution_gfname[nequation]);
          CCTK_REAL *elap = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, exact_laplacian_gfname[nequation]);

          CCTK_INT npoints = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];

          LC_LOOP3 (EL_CMP,
                    i, j, k, 0, 0, 0, cctk_lsh[0], cctk_lsh[1], cctk_lsh[2], 
                    cctk_lsh[0], cctk_lsh[1], cctk_lsh[2])
          {
            int sindex = CCTK_GFINDEX3D(cctkGH,i,j,k);
            int index = sindex + nequation * npoints;

            ct_terr[index] = epsi[sindex] - ct_psi[index] - exact_offset;
            if (outdated_psi) ct_terr[index] -= ct_err[index];
          } LC_ENDLOOP3 (EL_CMP);

          int xgh = cctk_nghostzones[0];
          int ygh = cctk_nghostzones[1];
          int zgh = cctk_nghostzones[2];
          int imin = xgh; 
          int imax = cctk_lsh[0]-xgh;
          int jmin = ygh;
          int jmax = cctk_lsh[1]-ygh;
          int kmin = zgh;
          int kmax = cctk_lsh[2]-zgh;

          LC_LOOP3 (EL_TRC,
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

            CCTK_REAL epsix, epsiy, epsiz, epsixx, epsiyy, epsizz, epsixy, epsixz, epsiyz;
            CT_FD(epsi, stnc, 0,
                  &epsix, &epsiy, &epsiz, &epsixx, &epsiyy, &epsizz, &epsixy, &epsixz, &epsiyz);

            ct_trunc[index] = elap[sindex] - epsixx - epsiyy - epsizz;
            ct_trunc_copy[index] = ct_trunc[index];
          } LC_ENDLOOP3 (EL_TRC);
        }

      } END_COMPONENT_LOOP;
    } END_MAP_LOOP;

    CCTK_REAL norm;
    CT_Norm(CCTK_PASS_CTOC, "CT_MultiLevel::ct_terr[0]", &norm, nequation);

    char fname[100];
    sprintf(fname, "terr_norm_eqn%d.asc", nequation);
    CT_WriteTimeSeries(Carpet::reflevel, norm, fname); 

    CT_Norm(CCTK_PASS_CTOC, "CT_MultiLevel::ct_trunc[0]", &norm, nequation);

    sprintf(fname, "trunc_norm_eqn%d.asc", nequation);
    CT_WriteTimeSeries(Carpet::reflevel, norm, fname); 
  }

  return;
}
