#include <fftw3.h>
#include <math.h>

int main(){
  int nx=2; int ny=5;
  int N = nx*ny;
  double *in;
  fftw_plan p;
  in = (double*) fftw_malloc(sizeof(double) * N);
  p = fftw_plan_r2r_2d(ny,nx,in,in,FFTW_RODFT00,FFTW_RODFT00,FFTW_MEASURE);

  double dx = 1./(nx+1); double xx = 0;
  double dy = 1./(ny+1); double yy = 0;
  int idx = 0;
  for(int j=0;j<ny;j++){
    for(int i=0;i<nx;i++){
      idx = j*nx+i;
      xx = i*dx;
      yy = j*dy;
      in[idx] = xx*xx + exp(yy);
    }
  }

  printf("Function Data\n");
  for(int i=0;i<N;i++) printf("%f\t",in[i]);
  printf("\n");

  fftw_execute(p);
  printf("Sine Transformed\n");
  for(int i=0;i<N;i++) printf("%f\t",in[i]);
  printf("\n");

  double d = 2;
  for(int i=0;i<N;i++) in[i] = (double) in[i]*d;
  printf("Did some manipulation\n");
  for(int i=0;i<N;i++) printf("%f\t",in[i]);
  printf("\n");

  fftw_execute(p); /* find correct factor */
  double factor = (double) 0.25/(ny+1)/(nx+1);
  for(int i=0;i<N;i++) in[i] = in[i]*factor;
  printf("Inverse Transformed\n");
  for(int i=0;i<N;i++) printf("%f\t",in[i]);
  printf("\n");

  fftw_destroy_plan(p);
  fftw_free(in);  
  return 0;
}
