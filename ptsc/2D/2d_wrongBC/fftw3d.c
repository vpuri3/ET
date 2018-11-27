#include <fftw3.h>
#include <math.h>

int main(){
  int nx=2, ny=2, nz=2;
  int N = nx*ny*nz;
  double d = 2;
  double *in;
  fftw_plan p;
  in = (double*) fftw_malloc(sizeof(double) * N);
  p = fftw_plan_r2r_3d(nz,ny,nx,in,in,FFTW_RODFT00,FFTW_RODFT00,FFTW_RODFT00,FFTW_MEASURE);
  double factor = (double) 0.125/(nz+1)/(ny+1)/(nx+1);

  double dx = 1./(nx+1); double xx = 0;
  double dy = 1./(ny+1); double yy = 0;
  double dz = 1./(nz+1); double zz = 0;
  int idx = 0;
  for(int k=0;k<nz;k++){
    for(int j=0;j<ny;j++){
      for(int i=0;i<nx;i++){
	idx = k*ny*nx + j*nx + i;
	xx = i*dx;
	yy = j*dy;
	zz = k*dz;
	in[idx] = xx+yy+zz;
      }
    }
  }

  printf("Function Data\n");
  for(int i=0;i<N;i++) printf("%f\t",in[i]);
  printf("\n");

  fftw_execute(p);
  printf("Sine Transformed\n");
  for(int i=0;i<N;i++) printf("%f\t",in[i]);
  printf("\n");

  for(int i=0;i<N;i++) in[i] = (double) in[i]*d;
  printf("Did some manipulation\n");
  for(int i=0;i<N;i++) printf("%f\t",in[i]);
  printf("\n");

  fftw_execute(p);
  for(int i=0;i<N;i++) in[i] = in[i]*factor;
  printf("Inverse Transformed\n");
  for(int i=0;i<N;i++) printf("%f\t",in[i]);
  printf("\n");

  fftw_destroy_plan(p);
  fftw_free(in);  
  return 0;
}
