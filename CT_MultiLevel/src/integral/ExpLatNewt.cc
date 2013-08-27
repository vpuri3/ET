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
CCTK_REAL A=0;

for (int n=1; n<10; n++)
{
  LC_LOOP3 (EL_ENF,
            i, j, k, 0, 0, 0, cctk_lsh[0], cctk_lsh[1], cctk_lsh[2], 
            cctk_lsh[0], cctk_lsh[1], cctk_lsh[2])
  {
    int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

    ct_integrand1[index] = ct_c2[index] * pow(ct_psi[index]+ct_a2[index]+A, n2[0]);
    ct_integrand2[index] = ct_c1[index] * pow(ct_psi[index]+ct_a1[index]+A, n1[0]);
    ct_integrand3[index] = n2[0] * ct_c2[index] * pow(ct_psi[index]+ct_a2[index]+A, n2[0]-1);
    ct_integrand4[index] = n1[0] * ct_c1[index] * pow(ct_psi[index]+ct_a1[index]+A, n1[0]-1);
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

  A -= ( int_1 + int_2 + 2*pi*mass ) / (int_3 + int_4);

//  cout << int_1 + int_2 + 2*pi*mass << "\t" << A << endl;
}

value[0] = A;

CT_WriteTimeSeries(Carpet::reflevel, A, "A.asc"); 

if (CCTK_Equals(veryverbose, "yes")) CCTK_VInfo(CCTK_THORNSTRING, "Evaluating the integral condition (reset value = %f)...", value[0]);
