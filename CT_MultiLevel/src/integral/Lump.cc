//FIXME: this below assumes that K is a constant
CCTK_REAL *K = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, "CT_Analytic::testK"); 
CCTK_REAL *c0 = (CCTK_REAL *) CCTK_VarDataPtr(cctkGH, 0, "CT_Analytic::testc0"); 

LC_LOOP3 (EL_ENF,
          i, j, k, 0, 0, 0, cctk_lsh[0], cctk_lsh[1], cctk_lsh[2], 
          cctk_lsh[0], cctk_lsh[1], cctk_lsh[2])
{
  int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

  ct_integrand1[index] = ct_c1[index] * pow(ct_psi[index]+ct_a1[index], n1[0]) + ct_c2[index] * pow(ct_psi[index]+ct_a2[index], n2[0]);
  ct_integrand2[index] = pow(ct_psi[index]+ct_a0[index], n0[0]) / 12.0;
} LC_ENDLOOP3 (EL_ENF);

CCTK_REAL int_1, int_2;
char name[100];
sprintf(name, "CT_MultiLevel::ct_integrand1[0]");
CT_Integral(CCTK_PASS_CTOC, name, &int_1, integral_refinement, 0);
sprintf(name, "CT_MultiLevel::ct_integrand2[0]");
CT_Integral(CCTK_PASS_CTOC, name, &int_2, integral_refinement, 0);

CCTK_REAL Kfin = -sqrt(abs(int_1/int_2)); // FIXME: The argument of abs shouldn't be negative. If it is, setting Kfin this way is wrong.
LC_LOOP3 (EL_ENF,
          i, j, k, 0, 0, 0, cctk_lsh[0], cctk_lsh[1], cctk_lsh[2], 
          cctk_lsh[0], cctk_lsh[1], cctk_lsh[2])
{
  int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

  K[index] = Kfin;
  c0[index] = - Kfin * Kfin / 12.0;

} LC_ENDLOOP3 (EL_ENF);

CT_WriteTimeSeries(Carpet::reflevel, Kfin, "k.asc"); 

if (CCTK_Equals(veryverbose, "yes")) CCTK_VInfo(CCTK_THORNSTRING, "Evaluating the integral condition (reset value = %f)...", value[0]);
