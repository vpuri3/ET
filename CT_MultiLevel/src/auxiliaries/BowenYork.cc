CCTK_INT nauxiliary;

nauxiliary = 0;
int aindex = sindex + nauxiliary * stnc.npoints;

CCTK_REAL *Axx = cptr.extras[0]; 
CCTK_REAL *Axy = cptr.extras[1];
CCTK_REAL *Axz = cptr.extras[2];
CCTK_REAL *Ayy = cptr.extras[3];
CCTK_REAL *Ayz = cptr.extras[4];
CCTK_REAL *Azz = cptr.extras[5];

CCTK_REAL L2 = Axx[sindex]*Axx[sindex] + 2*Axy[sindex]*Axy[sindex] + 2*Axz[sindex]*Axz[sindex] + Ayy[sindex]*Ayy[sindex] + 2*Ayz[sindex]*Ayz[sindex] + Azz[sindex]*Azz[sindex];

cptr.ct_auxiliary[aindex] = L2 / 8.0;
