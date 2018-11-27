#include <fftw3.h>
#include <math.h>

int main(){
  int N = 10;
  double *in;
  fftw_plan p;
  in = (double*) fftw_malloc(sizeof(double) * N);
  p = fftw_plan_r2r_1d(N,in,in,FFTW_RODFT00, FFTW_MEASURE);
  /* fill in array now */
  double dx = (double) 1/(N+1); double xx = 0;
  for(int i=0;i<N;i++){
    xx += dx;
    in[i] = xx*xx + exp(xx);
  }

  printf("Function Data\n");
  for(int i=0;i<N;i++) printf("%f\t",in[i]);
  printf("\n");

  fftw_execute(p);
  printf("Sine Transformed\n");
  for(int i=0;i<N;i++) printf("%f\t",in[i]);
  printf("\n");

  for(int i=0;i<N;i++) in[i]= 2*in[i];

  printf("Did some manipulation\n");
  for(int i=0;i<N;i++) printf("%f\t",in[i]);
  printf("\n");

  fftw_execute(p);
  double factor = (double) 0.5/(N+1);
  for(int i=0;i<N;i++) in[i] = in[i]*factor;
  printf("Inverse Transformed\n");
  for(int i=0;i<N;i++) printf("%f\t",in[i]);
  printf("\n");

  fftw_destroy_plan(p);
  fftw_free(in);  
  return 0;
}
