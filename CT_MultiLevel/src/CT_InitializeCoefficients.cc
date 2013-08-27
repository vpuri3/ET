#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "carpet.hh"
#include "loopcontrol.h"
#include "Carpet/Carpet/src/modes.hh"

#include "CT_MultiLevel.hh"

extern "C" void CT_InitializeCoefficients(int sindex, int nequation, struct coeffptr *cptr, struct stencil stnc)
{
  int index = sindex + nequation * stnc.npoints;

  CT_SetAuxiliaries(sindex, stnc, *cptr);

  cptr->ct_cxx[index] = cptr->cxx[sindex];
  cptr->ct_cxy[index] = cptr->cxy[sindex];
  cptr->ct_cxz[index] = cptr->cxz[sindex];
  cptr->ct_cyy[index] = cptr->cyy[sindex];
  cptr->ct_cyz[index] = cptr->cyz[sindex];
  cptr->ct_czz[index] = cptr->czz[sindex];
  cptr->ct_cx[index] = cptr->cx[sindex];
  cptr->ct_cy[index] = cptr->cy[sindex];
  cptr->ct_cz[index] = cptr->cz[sindex];
  cptr->ct_c0[index] = cptr->c0[sindex];
  cptr->ct_c1[index] = cptr->c1[sindex];
  cptr->ct_c2[index] = cptr->c2[sindex];
  cptr->ct_c3[index] = cptr->c3[sindex];
  cptr->ct_c4[index] = cptr->c4[sindex];
  cptr->ct_a0[index] = cptr->a0[sindex];
  cptr->ct_a1[index] = cptr->a1[sindex];
  cptr->ct_a2[index] = cptr->a2[sindex];
  cptr->ct_a3[index] = cptr->a3[sindex];
  cptr->ct_a4[index] = cptr->a4[sindex];

  return;
}
