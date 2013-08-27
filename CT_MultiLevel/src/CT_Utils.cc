#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "carpet.hh"
#include "loopcontrol.h"

#include "CT_MultiLevel.hh"

#ifdef CCTK_MPI
#  include <mpi.h>
#endif

using namespace Carpet;

int CT_ProcessOwnsData()
{
  int rank;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int pc = vhh.AT(Carpet::mglevel)->processor(Carpet::reflevel,Carpet::component);

  return (rank == pc);
}
  
extern "C" void CT_Restrict(CCTK_ARGUMENTS, const char *varname)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_Equals(verbose, "yes")) CCTK_VInfo(CCTK_THORNSTRING, "Restricting %s...", varname);

  int ml=Carpet::mglevel;
  int ll=Carpet::reflevel;
  for (comm_state state; not state.done(); state.step()) {
    const int g = CCTK_GroupIndex(varname); 
    const int active_tl = CCTK_ActiveTimeLevelsGI (cctkGH, g);
    assert (active_tl>=0);
    const int tl = active_tl > 1 ? timelevel : 0;
    for (int m=0; m<(int)arrdata.AT(g).size(); ++m) {
      for (int v = 0; v < (int)arrdata.AT(g).AT(m).data.size(); ++v) {
        ggf *const gv = arrdata.AT(g).AT(m).data.AT(v);
        gv->ref_restrict_all (state, tl, ll, ml);
      }
    }
  } // for state

  CT_UpdateBoundaries(CCTK_PASS_CTOC, varname);

  return;
}

extern "C" void CT_Prolongate(CCTK_ARGUMENTS, const char *varname)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_Equals(verbose, "yes")) CCTK_VInfo(CCTK_THORNSTRING, "Prolongating %s...", varname);

  int ml=Carpet::mglevel;
  int ll=Carpet::reflevel;
  for (comm_state state; not state.done(); state.step()) {
    const int g = CCTK_GroupIndex(varname); 
    const int active_tl = CCTK_ActiveTimeLevelsGI (cctkGH, g);
    assert (active_tl>=0);
    const int tl = active_tl > 1 ? timelevel : 0;
    for (int m=0; m<(int)arrdata.AT(g).size(); ++m) {
      for (int v = 0; v < (int)arrdata.AT(g).AT(m).data.size(); ++v) {
        ggf *const gv = arrdata.AT(g).AT(m).data.AT(v);
        gv->ref_prolongate_all (state, tl, ll, ml, 0);
        gv->ref_bnd_prolongate_all (state, tl, ll, ml, 0);
      }
    }
  } // for state

  CT_UpdateBoundaries(CCTK_PASS_CTOC, varname);

  return;
}

extern "C" void CT_ProlongateBndrs(CCTK_ARGUMENTS, const char *varname)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_Equals(verbose, "yes")) CCTK_VInfo(CCTK_THORNSTRING, "Prolongating boundaries of %s...", varname);

  int ml=Carpet::mglevel;
  int ll=Carpet::reflevel;
  for (comm_state state; not state.done(); state.step()) {
    const int g = CCTK_GroupIndex(varname); 
    const int active_tl = CCTK_ActiveTimeLevelsGI (cctkGH, g);
    assert (active_tl>=0);
    const int tl = active_tl > 1 ? timelevel : 0;
    for (int m=0; m<(int)arrdata.AT(g).size(); ++m) {
      for (int v = 0; v < (int)arrdata.AT(g).AT(m).data.size(); ++v) {
        ggf *const gv = arrdata.AT(g).AT(m).data.AT(v);
//        gv->ref_prolongate_all (state, tl, ll, ml, 0);
        gv->ref_bnd_prolongate_all (state, tl, ll, ml, 0);
      }
    }
  } // for state

  return;
}

extern "C" void CT_FD(CCTK_REAL *gfunc, struct stencil stnc, int nequation,
                      CCTK_REAL *derx, CCTK_REAL *dery, CCTK_REAL *derz,
                      CCTK_REAL *derxx, CCTK_REAL *deryy, CCTK_REAL *derzz,
                      CCTK_REAL *derxy, CCTK_REAL *derxz, CCTK_REAL *deryz)
{
  int shft = nequation * stnc.npoints;
  if (stnc.order == 2)
  {
    *derx  = ( gfunc[stnc.px+shft] - gfunc[stnc.mx+shft] ) / (2.0*stnc.dx);
    *dery  = ( gfunc[stnc.py+shft] - gfunc[stnc.my+shft] ) / (2.0*stnc.dy);
    *derz  = ( gfunc[stnc.pz+shft] - gfunc[stnc.mz+shft] ) / (2.0*stnc.dz);
    *derxx = ( gfunc[stnc.px+shft] + gfunc[stnc.mx+shft] - 2.0 * gfunc[stnc.c+shft] ) / (stnc.dx*stnc.dx);
    *deryy = ( gfunc[stnc.py+shft] + gfunc[stnc.my+shft] - 2.0 * gfunc[stnc.c+shft] ) / (stnc.dy*stnc.dy);
    *derzz = ( gfunc[stnc.pz+shft] + gfunc[stnc.mz+shft] - 2.0 * gfunc[stnc.c+shft] ) / (stnc.dz*stnc.dz);
    *derxy = ( gfunc[stnc.mxmy+shft] - gfunc[stnc.mxpy+shft] - gfunc[stnc.pxmy+shft] + gfunc[stnc.pxpy+shft] ) / (4*stnc.dx*stnc.dy);
    *derxz = ( gfunc[stnc.mxmz+shft] - gfunc[stnc.mxpz+shft] - gfunc[stnc.pxmz+shft] + gfunc[stnc.pxpz+shft] ) / (4*stnc.dx*stnc.dz);
    *deryz = ( gfunc[stnc.mymz+shft] - gfunc[stnc.mypz+shft] - gfunc[stnc.pymz+shft] + gfunc[stnc.pypz+shft] ) / (4*stnc.dy*stnc.dz);
  } else if (stnc.order == 4) {
    *derx  = ( gfunc[stnc.m2x+shft] - gfunc[stnc.p2x+shft] ) / (12.0*stnc.dx) + 2.0 * ( gfunc[stnc.px+shft] - gfunc[stnc.mx+shft] ) / (3.0*stnc.dx);
    *dery  = ( gfunc[stnc.m2y+shft] - gfunc[stnc.p2y+shft] ) / (12.0*stnc.dy) + 2.0 * ( gfunc[stnc.py+shft] - gfunc[stnc.my+shft] ) / (3.0*stnc.dy);
    *derz  = ( gfunc[stnc.m2z+shft] - gfunc[stnc.p2z+shft] ) / (12.0*stnc.dz) + 2.0 * ( gfunc[stnc.pz+shft] - gfunc[stnc.mz+shft] ) / (3.0*stnc.dz);
    *derxx = - ( gfunc[stnc.m2x+shft] + gfunc[stnc.p2x+shft] ) / (12.0*stnc.dx*stnc.dx) + 4.0 * (gfunc[stnc.px+shft] + gfunc[stnc.mx+shft]) / (3.0*stnc.dx*stnc.dx) - 2.5 * gfunc[stnc.c+shft] / (stnc.dx*stnc.dx);
    *deryy = - ( gfunc[stnc.m2y+shft] + gfunc[stnc.p2y+shft] ) / (12.0*stnc.dy*stnc.dy) + 4.0 * (gfunc[stnc.py+shft] + gfunc[stnc.my+shft]) / (3.0*stnc.dy*stnc.dy) - 2.5 * gfunc[stnc.c+shft] / (stnc.dy*stnc.dy);
    *derzz = - ( gfunc[stnc.m2z+shft] + gfunc[stnc.p2z+shft] ) / (12.0*stnc.dz*stnc.dz) + 4.0 * (gfunc[stnc.pz+shft] + gfunc[stnc.mz+shft]) / (3.0*stnc.dz*stnc.dz) - 2.5 * gfunc[stnc.c+shft] / (stnc.dz*stnc.dz);
    // FIXME: the expressions below are only 2nd order accurate!
    *derxy = ( gfunc[stnc.mxmy+shft] - gfunc[stnc.mxpy+shft] - gfunc[stnc.pxmy+shft] + gfunc[stnc.pxpy+shft] ) / (4*stnc.dx*stnc.dy);
    *derxz = ( gfunc[stnc.mxmz+shft] - gfunc[stnc.mxpz+shft] - gfunc[stnc.pxmz+shft] + gfunc[stnc.pxpz+shft] ) / (4*stnc.dx*stnc.dz);
    *deryz = ( gfunc[stnc.mymz+shft] - gfunc[stnc.mypz+shft] - gfunc[stnc.pymz+shft] + gfunc[stnc.pypz+shft] ) / (4*stnc.dy*stnc.dz);
  } else {
    CCTK_WARN(0, "Differencing order not supported");
  }

  return;
}

extern "C" void CT_PrintError(CCTK_ARGUMENTS, CCTK_INT nequation)
{
  BEGIN_MAP_LOOP(cctkGH,CCTK_GF) {
    BEGIN_COMPONENT_LOOP(cctkGH, CCTK_GF) {
      DECLARE_CCTK_ARGUMENTS;
      DECLARE_CCTK_PARAMETERS;

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

        LC_LOOP3 (EL_INI,
                  //i, j, k, 0, 0, 0, cctk_lsh[0], cctk_lsh[1], cctk_lsh[2], 
                  i, j, k, imin, jmin, kmin, imax, jmax, kmax,
                  cctk_lsh[0], cctk_lsh[1], cctk_lsh[2])
        {
          int index = CCTK_GFINDEX3D(cctkGH,i,j,k) + nequation * npoints;

          printf("%f\t", ct_err[index]);
        } LC_ENDLOOP3 (EL_INI);
      }

    } END_COMPONENT_LOOP;
  } END_MAP_LOOP;

  printf("\n");

  return;
}

extern "C" void CT_PrintPsi(CCTK_ARGUMENTS, CCTK_INT nequation)
{
  BEGIN_MAP_LOOP(cctkGH,CCTK_GF) {
    BEGIN_COMPONENT_LOOP(cctkGH, CCTK_GF) {
      DECLARE_CCTK_ARGUMENTS;
      DECLARE_CCTK_PARAMETERS;

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

        LC_LOOP3 (EL_INI,
                  //i, j, k, 0, 0, 0, cctk_lsh[0], cctk_lsh[1], cctk_lsh[2], 
                  i, j, k, imin, jmin, kmin, imax, jmax, kmax,
                  cctk_lsh[0], cctk_lsh[1], cctk_lsh[2])
        {
          int index = CCTK_GFINDEX3D(cctkGH,i,j,k) + nequation * npoints;

          printf("%f\t", ct_psi[index]);
        } LC_ENDLOOP3 (EL_INI);
      }

      CCTK_REAL norm;
      char psiname[100];
      sprintf(psiname, "CT_MultiLevel::ct_residual[0]");
      CT_Norm(CCTK_PASS_CTOC, psiname, &norm, nequation);
      printf("norm = %f\n", norm);
    } END_COMPONENT_LOOP;
  } END_MAP_LOOP;

  printf("\n");

  return;
}

//extern "C" void CT_Norm(CCTK_ARGUMENTS, const string *varName, CCTK_REAL *norm)
extern "C" void CT_Norm(CCTK_ARGUMENTS, char *varName, CCTK_REAL *norm, CCTK_INT nequation)
{
  BEGIN_MAP_LOOP(cctkGH,CCTK_GF) {
    BEGIN_COMPONENT_LOOP(cctkGH, CCTK_GF) {
      DECLARE_CCTK_ARGUMENTS;
      DECLARE_CCTK_PARAMETERS;

      if (CT_ProcessOwnsData())
      {
        CCTK_REAL *var = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, varName);

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

        int n = 0;
        CCTK_REAL sum = 0;

#pragma omp parallel reduction(+: sum, n)
        LC_LOOP3 (EL_NRM,
                  i, j, k, imin, jmin, kmin, imax, jmax, kmax,
                  cctk_lsh[0], cctk_lsh[1], cctk_lsh[2])
        {
          int index = CCTK_GFINDEX3D(cctkGH,i,j,k) + nequation * npoints;

          sum += fabs(var[index]);
          n++;
        } LC_ENDLOOP3 (EL_NRM);

        CCTK_REAL localnorm = sum / n;
        MPI_Reduce(&localnorm, norm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Bcast(norm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      } 
    } END_COMPONENT_LOOP;
  } END_MAP_LOOP;

  return;
}

extern "C" void CT_Integral(CCTK_ARGUMENTS, const char *varName, CCTK_REAL *integ, CCTK_INT refinement, CCTK_INT nequation)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (CT_ProcessOwnsData())
  {
    CCTK_REAL *var = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, varName);

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

    CCTK_REAL sum = 0, dV = CCTK_DELTA_SPACE(0) * CCTK_DELTA_SPACE(1) * CCTK_DELTA_SPACE(2);

#pragma omp parallel reduction(+: sum)
    LC_LOOP3 (EL_NRM,
              i, j, k, imin, jmin, kmin, imax, jmax, kmax,
              cctk_lsh[0], cctk_lsh[1], cctk_lsh[2])
    {
      int index = CCTK_GFINDEX3D(cctkGH,i,j,k) + nequation * npoints;

      sum += var[index] * dV;
    } LC_ENDLOOP3 (EL_NRM);

    CCTK_REAL localinteg = sum;
    MPI_Reduce(&localinteg, integ, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Bcast(integ, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  } 

  return;
}

/*extern "C" void CT_Integral(CCTK_ARGUMENTS, const char *varName, CCTK_REAL *integ, CCTK_INT reffac)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

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
    CCTK_REAL xmin = x[CCTK_GFINDEX3D(cctkGH,imin,jmin,kmin)];
    CCTK_REAL ymin = y[CCTK_GFINDEX3D(cctkGH,imin,jmin,kmin)];
    CCTK_REAL zmin = z[CCTK_GFINDEX3D(cctkGH,imin,jmin,kmin)];
    CCTK_REAL dx = CCTK_DELTA_SPACE(0)/(double)reffac, dy = CCTK_DELTA_SPACE(1)/(double)reffac, dz = CCTK_DELTA_SPACE(2)/(double)reffac;

    CCTK_INT npointsx = (imax-imin) * reffac; // Only for periodic domains
    CCTK_INT npointsy = (jmax-jmin) * reffac;
    CCTK_INT npointsz = (kmax-kmin) * reffac;
    CCTK_INT num_dims = 3, npoints=npointsx*npointsy*npointsz, num_input_arrays = 1, num_output_arrays = 1;
    CCTK_INT operator_handle = CCTK_InterpHandle("Hermite polynomial interpolation");
    CCTK_INT param_table_handle = Util_TableCreate(UTIL_TABLE_FLAGS_DEFAULT);
    CCTK_INT coord_system_handle = CCTK_CoordSystemHandle("cart3d");

    CCTK_REAL *xcoord = (CCTK_REAL*)malloc(npoints*sizeof(CCTK_REAL));
    CCTK_REAL *ycoord = (CCTK_REAL*)malloc(npoints*sizeof(CCTK_REAL));
    CCTK_REAL *zcoord = (CCTK_REAL*)malloc(npoints*sizeof(CCTK_REAL));

    for (int kdx=0; kdx<npointsz; kdx++)
      for (int jdx=0; jdx<npointsy; jdx++)
        for (int idx=0; idx<npointsx; idx++)
        {
          int index = kdx * (npointsx*npointsy) + jdx * npointsx + idx;
          xcoord[index] = xmin + idx * dx; 
          ycoord[index] = ymin + jdx * dy; 
          zcoord[index] = zmin + kdx * dz; 
        }

    CCTK_POINTER_TO_CONST interp_coords[3];
    interp_coords[0] = xcoord; 
    interp_coords[1] = ycoord; 
    interp_coords[2] = zcoord; 

    CCTK_INT input_array_indices[1] = {CCTK_VarIndex(varName)};
    CCTK_INT output_array_types[1] = {CCTK_VARIABLE_REAL};
    CCTK_REAL *output = (CCTK_REAL*)malloc(npoints*sizeof(CCTK_REAL));
    void *output_arrays[1] = {output};

    int ierr = Util_TableSetFromString(param_table_handle, "order=3");
    if (ierr < 0)
      CCTK_WARN(0, "Interpolation table couldn't be parsed.");

    ierr = CCTK_InterpGridArrays(cctkGH, num_dims, operator_handle, param_table_handle,
                                 coord_system_handle, npoints, CCTK_VARIABLE_REAL, interp_coords,
                                 num_input_arrays, input_array_indices, num_output_arrays,
                                 output_array_types, output_arrays);
    if (ierr < 0)
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING, "Interpolation of %s failed.", varName);

    CCTK_REAL sum = 0, dV = dx*dy*dz;

    for (int kdx=0; kdx<npointsz; kdx++)
      for (int jdx=0; jdx<npointsy; jdx++)
        for (int idx=0; idx<npointsx; idx++)
        {
          int index = kdx * (npointsx*npointsy) + jdx * npointsx + idx;

          sum += output[index] * dV;
        }

    CCTK_REAL localinteg = sum;
    MPI_Reduce(&localinteg, integ, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Bcast(integ, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    free(xcoord);
    free(ycoord);
    free(zcoord);
    free(output);
  } 

  return;
}*/

extern "C" void CT_FillFlagArrays(CCTK_INT *levels, CCTK_INT *downward, CCTK_INT toplevel)
{
  CCTK_INT down = 1;
  CCTK_INT nlevels = 2*(toplevel+1)-1;

  for (int irl=0; irl<nlevels; irl++)
  {
    levels[irl] = abs(toplevel-irl);

    if (down)
      downward[irl] = 1;
    else
      downward[irl] = 0;

    if (levels[irl] == 0) down = 0; 
  }

  return;
}

CCTK_REAL CT_GetValue(CCTK_ARGUMENTS, CCTK_REAL xcoord, CCTK_REAL ycoord, CCTK_REAL zcoord, const char *varName)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_REAL value;

  if (CT_ProcessOwnsData())
  {
    CCTK_INT num_dims = 3, npoints = 1, num_input_arrays = 1, num_output_arrays = 1; 
    CCTK_INT operator_handle = CCTK_InterpHandle("Hermite polynomial interpolation");
    CCTK_INT param_table_handle = Util_TableCreate(UTIL_TABLE_FLAGS_DEFAULT);
    CCTK_INT coord_system_handle = CCTK_CoordSystemHandle("cart3d");

    CCTK_POINTER_TO_CONST interp_coords[3];
    interp_coords[0] = &xcoord; 
    interp_coords[1] = &ycoord; 
    interp_coords[2] = &zcoord; 

    CCTK_INT input_array_indices[1] = {CCTK_VarIndex(varName)};
    CCTK_INT output_array_types[1] = {CCTK_VARIABLE_REAL};
    void *output_arrays[1] = { (void *) &value };

    int ierr = Util_TableSetFromString(param_table_handle, "order=3");
    if (ierr < 0)
      CCTK_WARN(0, "Interpolation table couldn't be parsed.");

    ierr = CCTK_InterpGridArrays(cctkGH, num_dims, operator_handle, param_table_handle,
                                 coord_system_handle, npoints, CCTK_VARIABLE_REAL, interp_coords,
                                 num_input_arrays, input_array_indices, num_output_arrays,
                                 output_array_types, output_arrays);
    if (ierr < 0)
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING, "Interpolation of %s failed.", varName);
  }

  return value;
}

extern "C" void CT_WriteTimeSeries(CCTK_INT iteration, CCTK_REAL value, char *filename)
{
  FILE *ofile;
  char fullname[500];
  const char *path;
  int rank;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0)
  {
    path = *((const char * const *) CCTK_ParameterGet("out_dir", "IOUtil", NULL));
    sprintf(fullname, "%s/%s", path, filename);

    ofile = fopen(fullname, "a");
    fprintf(ofile, "%d\t%1.12e\n", iteration, value);
    fclose(ofile);
  }

  return;
}

extern "C" void CT_WritePairs(CCTK_REAL value1, CCTK_REAL value2, char *filename)
{
  FILE *ofile;
  char fullname[500];
  const char *path;
  int rank;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0)
  {
    path = *((const char * const *) CCTK_ParameterGet("out_dir", "IOUtil", NULL));
    sprintf(fullname, "%s/%s", path, filename);

    ofile = fopen(fullname, "a");
    fprintf(ofile, "%1.12e\t%1.12e\n", value1, value2);
    fclose(ofile);
  }

  return;
}

extern "C" void CT_Copy(CCTK_ARGUMENTS, char *varName1, char *varName2, CCTK_INT nequation)
{
  BEGIN_MAP_LOOP(cctkGH,CCTK_GF) {
    BEGIN_COMPONENT_LOOP(cctkGH, CCTK_GF) {
      DECLARE_CCTK_ARGUMENTS;
      DECLARE_CCTK_PARAMETERS;

      if (CT_ProcessOwnsData())
      {
        CCTK_REAL *var1 = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, varName1);
        CCTK_REAL *var2 = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, varName2);

        CCTK_INT npoints = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];

        LC_LOOP3 (EL_INI,
                  i, j, k, 0, 0, 0, cctk_lsh[0], cctk_lsh[1], cctk_lsh[2], 
                  cctk_lsh[0], cctk_lsh[1], cctk_lsh[2])
        {
          int index = CCTK_GFINDEX3D(cctkGH,i,j,k) + nequation * npoints;

          var2[index] = var1[index];
        } LC_ENDLOOP3 (EL_INI);
      }

    } END_COMPONENT_LOOP;
  } END_MAP_LOOP;

  return;
}

extern "C" void CT_ClearFiles(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;

  char fullname[500];
  const char *path;

  path = *((const char * const *) CCTK_ParameterGet("out_dir", "IOUtil", NULL));

  for (int nequation=0; nequation < number_of_equations; nequation++)
  {
    sprintf(fullname, "%s/%s_eqn%d.asc", path, "psi_norm", nequation);
    remove(fullname);

    sprintf(fullname, "%s/%s_eqn%d.asc", path, "err_norm", nequation);
    remove(fullname);

    sprintf(fullname, "%s/%s_eqn%d.asc", path, "terr_norm", nequation);
    remove(fullname);

    sprintf(fullname, "%s/%s_eqn%d.asc", path, "trunc_norm", nequation);
    remove(fullname);

    sprintf(fullname, "%s/%s_eqn%d.asc", path, "walk", nequation);
    remove(fullname);
  }

  sprintf(fullname, "%s/k.asc", path);
  remove(fullname);

  sprintf(fullname, "%s/A.asc", path);
  remove(fullname);

  return;
}
