#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "carpet.hh"
#include "loopcontrol.h"
#include "Carpet/Carpet/src/modes.hh"

#include "CT_MultiLevel.hh"

extern "C" void FirstDerivative(CCTK_REAL *gfunc, struct stencil stnc, int nequation,
                                CCTK_REAL *psix, CCTK_REAL *psiy, CCTK_REAL *psiz);
extern "C" void SecondDerivative(CCTK_REAL *gfunc, struct stencil stnc, int nequation,
                                 CCTK_REAL *psixy, CCTK_REAL *psixz, CCTK_REAL *psiyz,
                                 CCTK_REAL *psixx, CCTK_REAL *psiyy, CCTK_REAL *psizz);
extern "C" void CT_CalculateA(struct stencil stnc,
                              CCTK_REAL *Xx, CCTK_REAL *Xy, CCTK_REAL *Xz, 
                              CCTK_REAL *Axx, CCTK_REAL *Axy, CCTK_REAL *Axz, 
                              CCTK_REAL *Ayy, CCTK_REAL *Ayz, CCTK_REAL *Azz, 
                              CCTK_REAL *A2);
extern "C" void CT_ConvertToBSSN(CCTK_ARGUMENTS);
extern "C" void CT_ConvertToCCZ4(CCTK_ARGUMENTS);

extern "C" void CT_SetAuxiliaries(int sindex, struct stencil stnc, struct coeffptr cptr)
{
  DECLARE_CCTK_PARAMETERS;

  if(CCTK_Equals(model, "Bowen-York"))
  {
    #include "auxiliaries/BowenYork.cc"
  }
  else if (CCTK_Equals(model, "Expanding BH lattice"))
  {
    #include "auxiliaries/ExpLat.cc"
  }
  else if (CCTK_Equals(model, "Lump"))
  {
    #include "auxiliaries/Lump.cc"
  }

  return;
}

extern "C" void CT_InitializeADM(CCTK_ARGUMENTS)
{
  BEGIN_REFLEVEL_LOOP(cctkGH) {
    BEGIN_MAP_LOOP(cctkGH,CCTK_GF) {
      BEGIN_COMPONENT_LOOP(cctkGH, CCTK_GF) {
        DECLARE_CCTK_ARGUMENTS;
        DECLARE_CCTK_PARAMETERS;

        if (CT_ProcessOwnsData())
        {
          if (CCTK_Equals(verbose, "yes")) CCTK_Info(CCTK_THORNSTRING, "Initializing the ADM variables...");

          for (int tl=0; tl<=2; tl++)
          {
            CCTK_REAL *admgxx = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ADMBase::gxx");
            CCTK_REAL *admgxy = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ADMBase::gxy");
            CCTK_REAL *admgxz = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ADMBase::gxz");
            CCTK_REAL *admgyy = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ADMBase::gyy");
            CCTK_REAL *admgyz = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ADMBase::gyz");
            CCTK_REAL *admgzz = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ADMBase::gzz");
            CCTK_REAL *admkxx = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ADMBase::kxx");
            CCTK_REAL *admkxy = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ADMBase::kxy");
            CCTK_REAL *admkxz = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ADMBase::kxz");
            CCTK_REAL *admkyy = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ADMBase::kyy");
            CCTK_REAL *admkyz = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ADMBase::kyz");
            CCTK_REAL *admkzz = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ADMBase::kzz");

            CCTK_REAL *trK = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, "CT_Analytic::testK");
            CCTK_REAL *a0  = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, "CT_Analytic::testa0");

            LC_LOOP3 (EL_INM,
                      i, j, k, 0, 0, 0, cctk_lsh[0], cctk_lsh[1], cctk_lsh[2], 
                      cctk_lsh[0], cctk_lsh[1], cctk_lsh[2])
            {
              int sindex = CCTK_GFINDEX3D(cctkGH,i,j,k);

              admgxx[sindex] = pow(ct_psi[sindex]+a0[sindex], 4.0);
              admgxy[sindex] = 0;
              admgxz[sindex] = 0;
              admgyy[sindex] = pow(ct_psi[sindex]+a0[sindex], 4.0);
              admgyz[sindex] = 0;
              admgzz[sindex] = pow(ct_psi[sindex]+a0[sindex], 4.0);
            } LC_ENDLOOP3 (EL_INM);
             
            if (CCTK_Equals(fill_Aij, "Analytic Aij"))
            {
	      CCTK_REAL *AxxE = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, "CT_Analytic::testAxx");
	      CCTK_REAL *AxyE = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, "CT_Analytic::testAxy");
	      CCTK_REAL *AxzE = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, "CT_Analytic::testAxz");
	      CCTK_REAL *AyyE = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, "CT_Analytic::testAyy");
	      CCTK_REAL *AyzE = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, "CT_Analytic::testAyz");
	      CCTK_REAL *AzzE = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, "CT_Analytic::testAzz");

              LC_LOOP3 (EL_INC,
                        i, j, k, 0, 0, 0, cctk_lsh[0], cctk_lsh[1], cctk_lsh[2], 
                        cctk_lsh[0], cctk_lsh[1], cctk_lsh[2])
              {
                int sindex = CCTK_GFINDEX3D(cctkGH,i,j,k);

	        CCTK_REAL Axx = AxxE[sindex];
	        CCTK_REAL Axy = AxyE[sindex];
	        CCTK_REAL Axz = AxzE[sindex];
	        CCTK_REAL Ayy = AyyE[sindex];
	        CCTK_REAL Ayz = AyzE[sindex];
	        CCTK_REAL Azz = AzzE[sindex];

  	        admkxx[sindex] = trK[sindex] * admgxx[sindex] / 3.0 + Axx * pow(ct_psi[sindex]+a0[sindex], -2.0);
  	        admkxy[sindex] = Axy * pow(ct_psi[sindex]+a0[sindex], -2.0);
	        admkxz[sindex] = Axz * pow(ct_psi[sindex]+a0[sindex], -2.0);
	        admkyy[sindex] = trK[sindex] * admgyy[sindex] / 3.0 + Ayy * pow(ct_psi[sindex]+a0[sindex], -2.0);
	        admkyz[sindex] = Ayz * pow(ct_psi[sindex]+a0[sindex], -2.0);
	        admkzz[sindex] = trK[sindex] * admgzz[sindex] / 3.0 + Azz * pow(ct_psi[sindex]+a0[sindex], -2.0);
              } LC_ENDLOOP3 (EL_INC);
            } else if (CCTK_Equals(fill_Aij, "Solver") || CCTK_Equals(fill_Aij, "Analytic Xi")) {
	      CCTK_REAL *Xx;
	      CCTK_REAL *Xy;
	      CCTK_REAL *Xz;

              if (CCTK_Equals(fill_Aij, "Solver"))
              {
	        Xx = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, "CT_MultiLevel::ct_psi[1]");
	        Xy = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, "CT_MultiLevel::ct_psi[2]");
	        Xz = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, "CT_MultiLevel::ct_psi[3]");
              } else if (CCTK_Equals(fill_Aij, "Analytic Xi")) {
	        Xx = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, "CT_Analytic::testXx");
	        Xy = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, "CT_Analytic::testXy");
	        Xz = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, "CT_Analytic::testXz");
              }

              LC_LOOP3 (EL_INC2,
                        i, j, k, 0, 0, 0, cctk_lsh[0], cctk_lsh[1], cctk_lsh[2], 
                        cctk_lsh[0], cctk_lsh[1], cctk_lsh[2])
              {
                int sindex = CCTK_GFINDEX3D(cctkGH,i,j,k);
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
                  cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]
                };

                CCTK_REAL Axx, Axy, Axz, Ayy, Ayz, Azz, A2;
                if (i>=cctk_nghostzones[0] && 
                    i<cctk_lsh[0]-cctk_nghostzones[0] && 
                    j>=cctk_nghostzones[1] && 
                    j<cctk_lsh[1]-cctk_nghostzones[1] && 
                    k>=cctk_nghostzones[2] && 
                    k<cctk_lsh[2]-cctk_nghostzones[2])
                {
                  CT_CalculateA(stnc, Xx, Xy, Xz, &Axx, &Axy, &Axz, &Ayy, &Ayz, &Azz, &A2);
                }
                else
                {
	          Axx = 0;
	          Axy = 0;
	          Axz = 0;
	          Ayy = 0;
	          Ayz = 0;
	          Azz = 0;
                }

  	        admkxx[sindex] = trK[sindex] * admgxx[sindex] / 3.0 + Axx * pow(ct_psi[sindex]+a0[sindex], -2.0);
  	        admkxy[sindex] = Axy * pow(ct_psi[sindex]+a0[sindex], -2.0);
	        admkxz[sindex] = Axz * pow(ct_psi[sindex]+a0[sindex], -2.0);
	        admkyy[sindex] = trK[sindex] * admgyy[sindex] / 3.0 + Ayy * pow(ct_psi[sindex]+a0[sindex], -2.0);
	        admkyz[sindex] = Ayz * pow(ct_psi[sindex]+a0[sindex], -2.0);
	        admkzz[sindex] = trK[sindex] * admgzz[sindex] / 3.0 + Azz * pow(ct_psi[sindex]+a0[sindex], -2.0);
              } LC_ENDLOOP3 (EL_INC2);
            } // if fill_Aij
          } // for tl
        } // if ProcessOwnsData
      } END_COMPONENT_LOOP;
    } END_MAP_LOOP;

    CT_UpdateBoundaries(CCTK_PASS_CTOC, "ADMBase::curv");

    if (CCTK_IsThornActive("ML_BSSN"))
    {
      BEGIN_MAP_LOOP(cctkGH,CCTK_GF) {
        BEGIN_COMPONENT_LOOP(cctkGH, CCTK_GF) {
        DECLARE_CCTK_ARGUMENTS;
        DECLARE_CCTK_PARAMETERS;

          if (CT_ProcessOwnsData())
          {
	    CT_ConvertToBSSN(CCTK_PASS_CTOC);
	    //CCTK_CallFunction("ML_BSSN_convertFromADMBase", NULL, NULL);
            CCTK_ScheduleTraverse("SetTmunu", cctkGH, NULL);
            CCTK_ScheduleTraverse("ML_BSSN_constraints1_group", cctkGH, NULL);
            CCTK_ScheduleTraverse("ML_BSSN_constraints2_group", cctkGH, NULL);
          }
        } END_COMPONENT_LOOP;
      } END_MAP_LOOP;
    }
    else if (CCTK_IsThornActive("ML_CCZ4"))
    {
      BEGIN_MAP_LOOP(cctkGH,CCTK_GF) {
        BEGIN_COMPONENT_LOOP(cctkGH, CCTK_GF) {
        DECLARE_CCTK_ARGUMENTS;
        DECLARE_CCTK_PARAMETERS;

          if (CT_ProcessOwnsData())
          {
	    CT_ConvertToCCZ4(CCTK_PASS_CTOC);
            CCTK_ScheduleTraverse("SetTmunu", cctkGH, NULL);
            CCTK_ScheduleTraverse("ML_CCZ4_constraints1_group", cctkGH, NULL);
            CCTK_ScheduleTraverse("ML_CCZ4_constraints2_group", cctkGH, NULL);
          }
        } END_COMPONENT_LOOP;
      } END_MAP_LOOP;
    }

  } END_REFLEVEL_LOOP;

  return;
}

extern "C" void FirstDerivative(CCTK_REAL *gfunc, struct stencil stnc, int nequation,
                                CCTK_REAL *psix, CCTK_REAL *psiy, CCTK_REAL *psiz)
{
  CCTK_REAL tmp;

  CT_FD(gfunc, stnc, nequation,
        psix, psiy, psiz, &tmp, &tmp, &tmp, &tmp, &tmp, &tmp);

  return;
}

extern "C" void SecondDerivative(CCTK_REAL *gfunc, struct stencil stnc, CCTK_INT nequation,
                                 CCTK_REAL *psixy, CCTK_REAL *psixz, CCTK_REAL *psiyz,
                                 CCTK_REAL *psixx, CCTK_REAL *psiyy, CCTK_REAL *psizz)
{
  CCTK_REAL tmp;

  CT_FD(gfunc, stnc, nequation,
        &tmp, &tmp, &tmp, psixx, psiyy, psizz, psixy, psixz, psiyz);

  return;
}

extern "C" void CT_CalculateA(struct stencil stnc,
                              CCTK_REAL *Xx, CCTK_REAL *Xy, CCTK_REAL *Xz, 
                              CCTK_REAL *Axx, CCTK_REAL *Axy, CCTK_REAL *Axz, 
                              CCTK_REAL *Ayy, CCTK_REAL *Ayz, CCTK_REAL *Azz, 
                              CCTK_REAL *A2)
{
  CCTK_REAL psiXx, psiXy, psiXz, psiYx, psiYy, psiYz, psiZx, psiZy, psiZz;

  FirstDerivative(Xx, stnc, 0, &psiXx, &psiXy, &psiXz);
  FirstDerivative(Xy, stnc, 0, &psiYx, &psiYy, &psiYz);
  FirstDerivative(Xz, stnc, 0, &psiZx, &psiZy, &psiZz);

  CCTK_REAL Z = psiXx + psiYy + psiZz;
  *Axx = 2.0 * (psiXx - Z / 3.0); 
  *Axy = (psiYx + psiXy);
  *Axz = (psiZx + psiXz);
  *Ayy = 2.0 * (psiYy - Z / 3.0); 
  *Ayz = (psiZy + psiYz);
  *Azz = 2.0 * (psiZz - Z / 3.0); 
  *A2 = (*Axx) * (*Axx) + 2.0 * (*Axy) * (*Axy) + 2.0 * (*Axz) * (*Axz) + (*Ayy) * (*Ayy) + 2.0 * (*Ayz) * (*Ayz) + (*Azz) * (*Azz);

  return;
}

extern "C" void CT_ConvertToBSSN(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_REAL one3rd = 1.0/3.0;

  for (int tl=0; tl<=2; tl++)
  {
    CCTK_REAL *bssna = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ML_BSSN::alpha");
    CCTK_REAL *bssnb1 = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ML_BSSN::beta1");
    CCTK_REAL *bssnb2 = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ML_BSSN::beta2");
    CCTK_REAL *bssnb3 = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ML_BSSN::beta3");

    CCTK_REAL *bssngxx = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ML_BSSN::gt11");
    CCTK_REAL *bssngxy = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ML_BSSN::gt12");
    CCTK_REAL *bssngxz = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ML_BSSN::gt13");
    CCTK_REAL *bssngyy = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ML_BSSN::gt22");
    CCTK_REAL *bssngyz = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ML_BSSN::gt23");
    CCTK_REAL *bssngzz = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ML_BSSN::gt33");

    CCTK_REAL *bssnphi = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ML_BSSN::phi");
    CCTK_REAL *bssntrK = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ML_BSSN::trK");

    CCTK_REAL *bssnAxx = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ML_BSSN::At11");
    CCTK_REAL *bssnAxy = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ML_BSSN::At12");
    CCTK_REAL *bssnAxz = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ML_BSSN::At13");
    CCTK_REAL *bssnAyy = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ML_BSSN::At22");
    CCTK_REAL *bssnAyz = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ML_BSSN::At23");
    CCTK_REAL *bssnAzz = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ML_BSSN::At33");

    CCTK_REAL *bssnXx = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ML_BSSN::Xt1");
    CCTK_REAL *bssnXy = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ML_BSSN::Xt2");
    CCTK_REAL *bssnXz = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ML_BSSN::Xt3");

    CCTK_REAL *admgxx = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ADMBase::gxx");
    CCTK_REAL *admgyy = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ADMBase::gyy");
    CCTK_REAL *admgzz = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ADMBase::gzz");
    CCTK_REAL *admkxx = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ADMBase::kxx");
    CCTK_REAL *admkxy = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ADMBase::kxy");
    CCTK_REAL *admkxz = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ADMBase::kxz");
    CCTK_REAL *admkyy = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ADMBase::kyy");
    CCTK_REAL *admkyz = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ADMBase::kyz");
    CCTK_REAL *admkzz = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ADMBase::kzz");

    // Assuming conformal flatness
    LC_LOOP3 (EL_CNV,
              i, j, k, 0, 0, 0, cctk_lsh[0], cctk_lsh[1], cctk_lsh[2], 
              cctk_lsh[0], cctk_lsh[1], cctk_lsh[2])
    {
      int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

      CCTK_REAL W = pow(admgxx[index], 0.25);
      CCTK_REAL traceK = W * W * (admkxx[index]+admkyy[index]+admkzz[index]); 

      bssna[index] = W;
      bssnb1[index] = 0;
      bssnb2[index] = 0;
      bssnb3[index] = 0;

      bssngxx[index] = 1;
      bssngxy[index] = 0;
      bssngxz[index] = 0;
      bssngyy[index] = 1;
      bssngyz[index] = 0;
      bssngzz[index] = 1;

      bssnphi[index] = W;
      bssntrK[index] = traceK;

      bssnAxx[index] = W * W * (admkxx[index] - one3rd * admgxx[index] * traceK);
      bssnAxy[index] = W * W * (admkxy[index]);
      bssnAxz[index] = W * W * (admkxz[index]);
      bssnAyy[index] = W * W * (admkyy[index] - one3rd * admgyy[index] * traceK);
      bssnAyz[index] = W * W * (admkyz[index]);
      bssnAzz[index] = W * W * (admkzz[index] - one3rd * admgzz[index] * traceK);

      bssnXx[index] = 0.0;
      bssnXy[index] = 0.0;
      bssnXz[index] = 0.0;

    } LC_ENDLOOP3 (EL_CNV);
  } // for tl

  return;
}

extern "C" void CT_ConvertToCCZ4(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_REAL one3rd = 1.0/3.0;

  for (int tl=0; tl<=2; tl++)
  {
    CCTK_REAL *ccz4a = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ML_CCZ4::alpha");
    CCTK_REAL *ccz4b1 = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ML_CCZ4::beta1");
    CCTK_REAL *ccz4b2 = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ML_CCZ4::beta2");
    CCTK_REAL *ccz4b3 = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ML_CCZ4::beta3");

    CCTK_REAL *ccz4gxx = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ML_CCZ4::gt11");
    CCTK_REAL *ccz4gxy = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ML_CCZ4::gt12");
    CCTK_REAL *ccz4gxz = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ML_CCZ4::gt13");
    CCTK_REAL *ccz4gyy = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ML_CCZ4::gt22");
    CCTK_REAL *ccz4gyz = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ML_CCZ4::gt23");
    CCTK_REAL *ccz4gzz = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ML_CCZ4::gt33");

    CCTK_REAL *ccz4phi = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ML_CCZ4::phi");
    CCTK_REAL *ccz4trK = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ML_CCZ4::trK");

    CCTK_REAL *ccz4Axx = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ML_CCZ4::At11");
    CCTK_REAL *ccz4Axy = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ML_CCZ4::At12");
    CCTK_REAL *ccz4Axz = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ML_CCZ4::At13");
    CCTK_REAL *ccz4Ayy = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ML_CCZ4::At22");
    CCTK_REAL *ccz4Ayz = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ML_CCZ4::At23");
    CCTK_REAL *ccz4Azz = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ML_CCZ4::At33");

    CCTK_REAL *ccz4Xx = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ML_CCZ4::Xt1");
    CCTK_REAL *ccz4Xy = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ML_CCZ4::Xt2");
    CCTK_REAL *ccz4Xz = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ML_CCZ4::Xt3");

    CCTK_REAL *admgxx = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ADMBase::gxx");
    CCTK_REAL *admgyy = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ADMBase::gyy");
    CCTK_REAL *admgzz = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ADMBase::gzz");
    CCTK_REAL *admkxx = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ADMBase::kxx");
    CCTK_REAL *admkxy = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ADMBase::kxy");
    CCTK_REAL *admkxz = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ADMBase::kxz");
    CCTK_REAL *admkyy = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ADMBase::kyy");
    CCTK_REAL *admkyz = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ADMBase::kyz");
    CCTK_REAL *admkzz = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, tl, "ADMBase::kzz");

    // Assuming conformal flatness
    LC_LOOP3 (EL_CNV,
              i, j, k, 0, 0, 0, cctk_lsh[0], cctk_lsh[1], cctk_lsh[2], 
              cctk_lsh[0], cctk_lsh[1], cctk_lsh[2])
    {
      int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

      CCTK_REAL W = pow(admgxx[index], 0.25);
      CCTK_REAL traceK = W * W * (admkxx[index]+admkyy[index]+admkzz[index]); 

      ccz4a[index] = W;
      ccz4b1[index] = 0;
      ccz4b2[index] = 0;
      ccz4b3[index] = 0;

      ccz4gxx[index] = 1;
      ccz4gxy[index] = 0;
      ccz4gxz[index] = 0;
      ccz4gyy[index] = 1;
      ccz4gyz[index] = 0;
      ccz4gzz[index] = 1;

      ccz4phi[index] = W;
      ccz4trK[index] = traceK;

      ccz4Axx[index] = W * W * (admkxx[index] - one3rd * admgxx[index] * traceK);
      ccz4Axy[index] = W * W * (admkxy[index]);
      ccz4Axz[index] = W * W * (admkxz[index]);
      ccz4Ayy[index] = W * W * (admkyy[index] - one3rd * admgyy[index] * traceK);
      ccz4Ayz[index] = W * W * (admkyz[index]);
      ccz4Azz[index] = W * W * (admkzz[index] - one3rd * admgzz[index] * traceK);

      ccz4Xx[index] = 0.0;
      ccz4Xy[index] = 0.0;
      ccz4Xz[index] = 0.0;

    } LC_ENDLOOP3 (EL_CNV);
  } // for tl

  return;
}
