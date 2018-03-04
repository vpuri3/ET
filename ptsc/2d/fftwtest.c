#include <fftw3.h>
#include <stdlib.h>
#include <math.h>

int main(){
  int N = 10;
  double *in, *out;
  fftw_plan p;
  in = (double*) fftw_malloc(sizeof(double) * N);
  out = (double*) fftw_malloc(sizeof(double) * N);
  //p = fftw_plan_dft_1d(N, in, out, FFTW_DFT, FFTW_MEASURE);
  p = fftw_plan_r2r_1d(N,in,out,FFTW_RODFT00, FFTW_MEASURE);
  /* fill in array now */
  fftw_execute(p);
  fftw_destroy_plan(p);
  fftw_free(in); fftw_free(out);
  
  return 0;
}
