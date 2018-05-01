//FIXME: this below assumes that K is a constant
CCTK_REAL *K = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, "CT_Analytic::testK"); 
CCTK_REAL *dxK = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, "CT_Analytic::testdxK"); 
CCTK_REAL *dyK = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, "CT_Analytic::testdyK"); 
CCTK_REAL *dzK = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, "CT_Analytic::testdzK"); 
CCTK_REAL *c1 = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, "CT_Analytic::testc1"); 
CCTK_REAL *c2 = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, "CT_Analytic::testc2"); 
CCTK_REAL *W = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, "CT_Analytic::testW"); 
CCTK_REAL *Kfac = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, "CT_Analytic::testc2"); 

CCTK_REAL Kini = Kfac[CCTK_GFINDEX3D(cctkGH,0,0,0)];
CCTK_REAL *massa = (CCTK_REAL *) CCTK_ParameterGet("massa", "CT_Analytic", NULL); 
CCTK_REAL *massb = (CCTK_REAL *) CCTK_ParameterGet("massb", "CT_Analytic", NULL); 
CCTK_REAL mass = *massa + *massb;
CCTK_REAL pi = acos(-1);

LC_LOOP3 (EL_ENF,
          i, j, k, 0, 0, 0, cctk_lsh[0], cctk_lsh[1], cctk_lsh[2], 
          cctk_lsh[0], cctk_lsh[1], cctk_lsh[2])
{
  int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

  ct_integrand1[index] = ct_c2[index] * pow(ct_psi[index]+ct_a2[index], n2[0]);
  ct_integrand2[index] = pow(W[index], 2.0) * pow(ct_psi[index]+ct_a1[index], n1[0]) / 12.0;
} LC_ENDLOOP3 (EL_ENF);

CCTK_REAL int_1, int_2;
char name[100];
sprintf(name, "CT_MultiLevel::ct_integrand1[0]");
CT_Integral(CCTK_PASS_CTOC, name, &int_1, integral_refinement, 0);
sprintf(name, "CT_MultiLevel::ct_integrand2[0]");
CT_Integral(CCTK_PASS_CTOC, name, &int_2, integral_refinement, 0);

CCTK_REAL Kfin = -sqrt(abs((int_1+2*pi*mass)/int_2)); // FIXME: The argument of abs shouldn't be negative. If it is, setting Kfin this way is wrong.
CCTK_REAL scale = Kfin / Kini;
LC_LOOP3 (EL_ENF,
          i, j, k, 0, 0, 0, cctk_lsh[0], cctk_lsh[1], cctk_lsh[2], 
          cctk_lsh[0], cctk_lsh[1], cctk_lsh[2])
{
  int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

  K[index] *= scale;
  dxK[index] *= scale;
  dyK[index] *= scale;
  dzK[index] *= scale;
  c1[index] *= scale*scale;
  c2[index] *= scale;

} LC_ENDLOOP3 (EL_ENF);

CT_WriteTimeSeries(Carpet::reflevel, Kfin, "k.asc"); 

if (CCTK_Equals(veryverbose, "yes")) CCTK_VInfo(CCTK_THORNSTRING, "Evaluating the integral condition (reset value = %f)...", value[0]);
