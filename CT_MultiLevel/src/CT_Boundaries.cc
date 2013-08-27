#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include "cctk_Loop.h"

#include "carpet.hh"
#include "Carpet/Carpet/src/modes.hh"

#include "CT_MultiLevel.hh"

extern "C" void CT_Boundaries(CCTK_ARGUMENTS, const char *varname)
{
  BEGIN_MAP_LOOP(cctkGH,CCTK_GF) {
    BEGIN_COMPONENT_LOOP(cctkGH, CCTK_GF) {
      DECLARE_CCTK_ARGUMENTS;
      DECLARE_CCTK_PARAMETERS;

      if (CT_ProcessOwnsData())
      {
        if (CCTK_Equals(verbose, "yes")) CCTK_Info(CCTK_THORNSTRING, "Apply BCs...");

        CCTK_INT npoints = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];

        int ivari = CCTK_FirstVarIndex(varname);
        int ivarf = CCTK_FirstVarIndex(varname) + CCTK_NumVarsInGroup(varname);

        CCTK_REAL c1 = - 0.5; 
        CCTK_REAL c2 =   2.0; 
        CCTK_REAL c3 = - 1.5; 
        CCTK_REAL idx = 1 / CCTK_DELTA_SPACE(0); 
        CCTK_REAL idy = 1 / CCTK_DELTA_SPACE(1); 
        CCTK_REAL idz = 1 / CCTK_DELTA_SPACE(2); 

        CCTK_INT bndsize[6], is_ghostbnd[6], is_symbnd[6], is_physbnd[6];
        GetBoundarySizesAndTypes(cctkGH, 6, bndsize, is_ghostbnd, is_symbnd, is_physbnd);
        CCTK_INT ivertex = CCTK_GFINDEX3D(cctkGH,cctk_lsh[0]-bndsize[1]-1,cctk_lsh[1]-bndsize[3]-1,cctk_lsh[2]-bndsize[5]-1);

        //Insert a check somewhere that TwoPunctures is active and puncture_u has storage
        CCTK_REAL *tpu = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, "TwoPunctures::puncture_u");

        for (int ivar=ivari; ivar<ivarf; ivar++)
        {
          CCTK_REAL *pvar = (CCTK_REAL *) CCTK_VarDataPtrI(cctkGH, 0, ivar);
          CCTK_REAL coeff = (pvar[ivertex] - 1) * r[ivertex];

          if (CCTK_Equals(boundary_conditions, "Robin"))
          {
            CCTK_LOOP3_BND(robin, cctkGH,
                           i, j, k, ni, nj, nk)
            {
              int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
  
              CCTK_REAL dirx = cctki2_idir;
              CCTK_REAL diry = cctki2_jdir;
              CCTK_REAL dirz = cctki2_kdir;
 
              int i1 = CCTK_GFINDEX3D(cctkGH,i-dirx,j,k);
              int j1 = CCTK_GFINDEX3D(cctkGH,i,j-diry,k);
              int k1 = CCTK_GFINDEX3D(cctkGH,i,j,k-dirz);
              int i2 = CCTK_GFINDEX3D(cctkGH,i-2*dirx,j,k);
              int j2 = CCTK_GFINDEX3D(cctkGH,i,j-2*diry,k);
              int k2 = CCTK_GFINDEX3D(cctkGH,i,j,k-2*dirz);
              int i3 = CCTK_GFINDEX3D(cctkGH,i-3*dirx,j,k);
              int j3 = CCTK_GFINDEX3D(cctkGH,i,j-3*diry,k);
              int k3 = CCTK_GFINDEX3D(cctkGH,i,j,k-3*dirz);
 
              CCTK_REAL pvarx = idx * (c1 * pvar[i1] + c2 * pvar[i2] + c3 * pvar[i3]); 
              CCTK_REAL pvary = idy * (c1 * pvar[j1] + c2 * pvar[j2] + c3 * pvar[j3]); 
              CCTK_REAL pvarz = idz * (c1 * pvar[k1] + c2 * pvar[k2] + c3 * pvar[k3]); 
 
              pvar[index] = 1 + coeff / r[index];//- /*r[index] */ ( dirx * pvarx + diry * pvary + dirz * pvarz);
 
            } CCTK_ENDLOOP3_BND (robin);
          } else if (CCTK_Equals(boundary_conditions, "TwoPunctures")) {
            //Insert a check somewhere that TwoPunctures is active and puncture_u has storage
            CCTK_REAL *tpu = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, "TwoPunctures::puncture_u");

            CCTK_LOOP3_BND(tp, cctkGH,
                           i, j, k, ni, nj, nk)
            {
              int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

              CCTK_REAL dirx = cctki2_idir;
              CCTK_REAL diry = cctki2_jdir;
              CCTK_REAL dirz = cctki2_kdir;

              int i1 = CCTK_GFINDEX3D(cctkGH,i-dirx,j,k);
              int j1 = CCTK_GFINDEX3D(cctkGH,i,j-diry,k);
              int k1 = CCTK_GFINDEX3D(cctkGH,i,j,k-dirz);
              int i2 = CCTK_GFINDEX3D(cctkGH,i-2*dirx,j,k);
              int j2 = CCTK_GFINDEX3D(cctkGH,i,j-2*diry,k);
              int k2 = CCTK_GFINDEX3D(cctkGH,i,j,k-2*dirz);
              int i3 = CCTK_GFINDEX3D(cctkGH,i-3*dirx,j,k);
              int j3 = CCTK_GFINDEX3D(cctkGH,i,j-3*diry,k);
              int k3 = CCTK_GFINDEX3D(cctkGH,i,j,k-3*dirz);

              CCTK_REAL pvarx = idx * (c1 * pvar[i1] + c2 * pvar[i2] + c3 * pvar[i3]); 
              CCTK_REAL pvary = idy * (c1 * pvar[j1] + c2 * pvar[j2] + c3 * pvar[j3]); 
              CCTK_REAL pvarz = idz * (c1 * pvar[k1] + c2 * pvar[k2] + c3 * pvar[k3]); 

              pvar[index] = 1.0 + tpu[index];

            } CCTK_ENDLOOP3_BND (tp);
          }
        }
      }

    } END_COMPONENT_LOOP;
  } END_MAP_LOOP;

  return;
}
