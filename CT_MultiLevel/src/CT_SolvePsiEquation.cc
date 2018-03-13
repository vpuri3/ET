#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "carpet.hh"
#include "loopcontrol.h"
#include "Carpet/Carpet/src/modes.hh"

#include "CT_MultiLevel.hh"

// this needs PETSc in ther REQUIRES line of configuration.ccl
// #include "petsc.h"

static void CT_PETScSolve(CCTK_ARGUMENTS);

extern "C" void CT_SolvePsiEquation(CCTK_ARGUMENTS, CCTK_REAL *norm, CCTK_INT *step)
{
  DECLARE_CCTK_PARAMETERS;

  BEGIN_MAP_LOOP(cctkGH,CCTK_GF) {
    BEGIN_COMPONENT_LOOP(cctkGH, CCTK_GF) {
      DECLARE_CCTK_ARGUMENTS;

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
        CCTK_REAL dtime = 0.5 / (pow(CCTK_DELTA_SPACE(0),-2.0) + pow(CCTK_DELTA_SPACE(1),-2.0) + pow(CCTK_DELTA_SPACE(2),-2.0)); 
        CCTK_REAL twothirds = 2.0/3.0;

        struct coeffptr *cptr = (struct coeffptr *) malloc(number_of_equations*sizeof(struct coeffptr));
        CT_PopulatePointerStruct(CCTK_PASS_CTOC, cptr);

        LC_LOOP3 (EL_RLX,
                  i, j, k, imin, jmin, kmin, imax, jmax, kmax,
                  cctk_lsh[0], cctk_lsh[1], cctk_lsh[2])
        {
          for (int nequation=0; nequation < number_of_equations; nequation++)
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
            // RH: the code above makes copies of the coefficients in
            // CT_MultiLevel internal variables. No idea what they are used for
            // but let's keep them for now.

#if 0
            CCTK_REAL psix, psiy, psiz, psixx, psiyy, psizz, psixy, psixz, psiyz;
            CT_FD(ct_psi, stnc, nequation, 
                  &psix, &psiy, &psiz, &psixx, &psiyy, &psizz, &psixy, &psixz, &psiyz);

            CCTK_REAL psi0 = pow(ct_psi[index]+ct_a0[index], n0[nequation]);
            CCTK_REAL psi1 = pow(ct_psi[index]+ct_a1[index], n1[nequation]);
            CCTK_REAL psi2 = pow(ct_psi[index]+ct_a2[index], n2[nequation]);
            CCTK_REAL psi3 = pow(ct_psi[index]+ct_a3[index], n3[nequation]);
            CCTK_REAL psi4 = pow(ct_psi[index]+ct_a4[index], n4[nequation]);

/*            ct_rhs[index]  =                 ( ct_cxx[index] * psixx
                                             + ct_cyy[index] * psiyy
                                             + ct_czz[index] * psizz
                                             + ct_cxy[index] * psixy
                                             + ct_cxz[index] * psixz
                                             + ct_cyz[index] * psiyz
                                             + ct_cx[index] * psix
                                             + ct_cy[index] * psiy 
                                             + ct_cz[index] * psiz
                                             + ct_c0[index] * psi0
                                             + ct_c1[index] * psi1
                                             + ct_c2[index] * psi2
                                             + ct_c3[index] * psi3
                                             + ct_c4[index] * psi4 );*/
            ct_psi[index] += omega * dtime * ( ct_cxx[index] * psixx
                                             + ct_cyy[index] * psiyy
                                             + ct_czz[index] * psizz
                                             + ct_cxy[index] * psixy
                                             + ct_cxz[index] * psixz
                                             + ct_cyz[index] * psiyz
                                             + ct_cx[index] * psix
                                             + ct_cy[index] * psiy 
                                             + ct_cz[index] * psiz
                                             + ct_c0[index] * psi0
                                             + ct_c1[index] * psi1
                                             + ct_c2[index] * psi2
                                             + ct_c3[index] * psi3
                                             + ct_c4[index] * psi4 
                                             );
#endif
          }
        } LC_ENDLOOP3 (EL_RLX);

        CT_PETScSolve(CCTK_PASS_CTOC);

        free(cptr);
      } // if
    } END_COMPONENT_LOOP;
  } END_MAP_LOOP;

  char name[100];
  *norm = 0;
  for (int nequation=0; nequation < number_of_equations; nequation++)
  {
    CCTK_REAL enorm;
    sprintf(name, "CT_MultiLevel::ct_psi[0]");
    CT_Norm(CCTK_PASS_CTOC, name, &enorm, nequation);

    if (nequation == 0) *step = *step + 1;

    if (CCTK_Equals(veryverbose,"yes")) CCTK_VInfo(CCTK_THORNSTRING, "Iteration %d, psinorm = %1.10e", *step, enorm);

    sprintf(name, "psi_norm_eqn%d.asc", nequation);
    if (CCTK_Equals(output_norms,"yes")) CT_WriteTimeSeries(*step, enorm, name); 
    *norm += enorm;
  }

  return;
}

static void CT_PETScSolve(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS; // gets me everything in param.ccl
  DECLARE_CCTK_ARGUMENTS; // gets me everything in interface.ccl

  int xgh = cctk_nghostzones[0];
  int ygh = cctk_nghostzones[1];
  int zgh = cctk_nghostzones[2];
  int imin = xgh;
  int imax = cctk_lsh[0]-xgh;
  int jmin = ygh;
  int jmax = cctk_lsh[1]-ygh;
  int kmin = zgh;
  int kmax = cctk_lsh[2]-zgh;

  // output some debug messages so that we know what is going on. Needs to be
  // to stderr so that they all show up.
  int myproc =  CCTK_MyProc(cctkGH);
  fprintf(stderr, "Solving on reflevel %d component %d iteration %d on process %d\n"
                  "Showing myproc, i,j,k, x,y,z, nequation, psi, ct_c0\n",
          Carpet::reflevel, Carpet::component, cctk_iteration, myproc);

  // number of points on the grid is nx * ny * nz
  int npoints = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];

  // RH: just a demo code to show how to access things. Your code must provide
  // updated values for the unknown field ct_psi in the range
  // xgh <= i < cctk_lsh[0] - xgh etc.
  for(int nequation = 0 ; nequation < number_of_equations ; nequation++) {
    for(int k = kmin ; k < kmax ; k++) {
      for(int j = jmin ; j < jmax ; j++) {
        for(int i = imin ; i < imax ; i++) {
          // linear index into a 3d array
          int index3d = i + j*cctk_lsh[0] + k*cctk_lsh[0]*cctk_lsh[1];
          // linear index into an vector of arrays eg psi(i,j,k,nequation) in matlab
          // that is the coefficients for the n'th equation (coutning from 0)
          // start at linear index: nequation * npoints
          int index4d = i + j*cctk_lsh[0] + k*cctk_lsh[0]*cctk_lsh[1] +
                        nequation*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2];

          // you cannot use fprintf(stdout, ...) for output since only one
          // process will output, but stderr will work.
          fprintf(stderr, "%d %d %d %d %g %g %g %d %g %g\n",
                  myproc, i,j,k, x[index3d],y[index3d],z[index3d], nequation,
                  ct_psi[index4d], ct_c0[index4d]);
        }
      }
    }
  }

  //int foo;
  //PetscOptionsGetInt(NULL,NULL,"-M",&foo,NULL);

  fprintf(stderr, "Done on reflevel %d component %d iteration %d on process %d\n",
          Carpet::reflevel, Carpet::component, cctk_iteration, myproc);
}
