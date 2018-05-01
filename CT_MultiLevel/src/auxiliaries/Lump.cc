int nauxiliary = 0;
int aindex = sindex + nauxiliary * stnc.npoints;

CCTK_REAL *Xx = &cptr.ct_psi[stnc.npoints]; 
CCTK_REAL *Xy = &cptr.ct_psi[2*stnc.npoints];
CCTK_REAL *Xz = &cptr.ct_psi[3*stnc.npoints];
CCTK_REAL *dxK = cptr.extras[0];
CCTK_REAL *dyK = cptr.extras[1];
CCTK_REAL *dzK = cptr.extras[2];

CCTK_REAL Lxx, Lxy, Lxz, Lyy, Lyz, Lzz, L2;
CT_CalculateA(stnc, Xx, Xy, Xz, &Lxx, &Lxy, &Lxz, &Lyy, &Lyz, &Lzz, &L2);

cptr.ct_auxiliary[aindex] = L2 / 8.0;

//psi^6 dK term in momentum constraint
CCTK_REAL psiXxy, psiXxz, psiYxy, psiYyz, psiZxz, psiZyz, psiXxx, psiYyy, psiZzz, tmp;
SecondDerivative(Xx, stnc, 0, &psiXxy, &psiXxz, &tmp, &psiXxx, &tmp, &tmp);
SecondDerivative(Xy, stnc, 0, &psiYxy, &tmp, &psiYyz, &tmp, &psiYyy, &tmp);
SecondDerivative(Xz, stnc, 0, &tmp, &psiZxz, &psiZyz, &tmp, &tmp, &psiZzz);

nauxiliary = 1;
aindex = sindex + nauxiliary * stnc.npoints;
cptr.ct_auxiliary[aindex] = ( -2.0 * pow(cptr.ct_psi[sindex], 6.0) * dxK[sindex] + psiXxx + psiYxy + psiZxz) / 3.0;

nauxiliary = 2;
aindex = sindex + nauxiliary * stnc.npoints;
cptr.ct_auxiliary[aindex] = ( -2.0 * pow(cptr.ct_psi[sindex], 6.0) * dyK[sindex] + psiXxy + psiYyy + psiZyz) / 3.0;

nauxiliary = 3;
aindex = sindex + nauxiliary * stnc.npoints;
cptr.ct_auxiliary[aindex] = ( -2.0 * pow(cptr.ct_psi[sindex], 6.0) * dzK[sindex] + psiXxz + psiYyz + psiZzz) / 3.0;
