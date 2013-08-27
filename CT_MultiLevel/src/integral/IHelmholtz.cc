LC_LOOP3 (EL_ENF,
          i, j, k, 0, 0, 0, cctk_lsh[0], cctk_lsh[1], cctk_lsh[2], 
          cctk_lsh[0], cctk_lsh[1], cctk_lsh[2])
{
  int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

  ct_integrand1[index] = ct_c0[index];
  ct_integrand2[index] = ct_c1[index];
  ct_integrand3[index] = ct_psi[index];
  ct_integrand4[index] = 1;
} LC_ENDLOOP3 (EL_ENF);

CCTK_REAL int_1, int_2, int_3, int_4;
char name[100];
sprintf(name, "CT_MultiLevel::ct_integrand1[0]");
CT_Integral(CCTK_PASS_CTOC, name, &int_1, integral_refinement, 0);
sprintf(name, "CT_MultiLevel::ct_integrand2[0]");
CT_Integral(CCTK_PASS_CTOC, name, &int_2, integral_refinement, 0);
sprintf(name, "CT_MultiLevel::ct_integrand3[0]");
CT_Integral(CCTK_PASS_CTOC, name, &int_3, integral_refinement, 0);
sprintf(name, "CT_MultiLevel::ct_integrand4[0]");
CT_Integral(CCTK_PASS_CTOC, name, &int_4, integral_refinement, 0);

value[0] = -int_1/int_2 - int_3/int_4;

if (CCTK_Equals(veryverbose, "yes")) CCTK_VInfo(CCTK_THORNSTRING, "Evaluating the integral condition (reset value = %f)...", value[0]);
