static char help[] = "Solves a tridiagonal linear system with KSP.\n\n";

#include <petscksp.h>
#include <math.h>
#include <fftw3.h>
#include "2d.h"

int main(int argc,char **args)
{
  PetscErrorCode ierr;
  Vec         x, u, b, d;
  Mat         A, F; /* Second order Laplace Operator, fft matrix */
  PetscInt    idx, nx=2, ny=5, N=nx*ny;
  Ctx         ctx;
  fftw_plan   p;

  ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,NULL,"-nx",&nx,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-ny",&ny,NULL);CHKERRQ(ierr);

  PetscReal * u_ = malloc(N*sizeof(PetscReal));
  PetscReal * x_ = malloc(N*sizeof(PetscReal));
  PetscReal * d_ = malloc(N*sizeof(PetscReal));

  ctx.nx = nx; ctx.ny = ny;
  ctx.dxinv = (double) nx+1; ctx.dyinv = (double) ny+1;
  ctx.dx = 1/ctx.dxinv; ctx.dy = 1/ctx.dyinv;
  ctx.p = p; ctx.fft_factor = (double) 0.25/(N+1)/(nx+1);
  ctx.x_ = x_;  
  ctx.d = &d;
  /*----------------CREATING VECTORS--------------*/
  for(int i=0;i<N;i++){ u_[i]=i; x_[i]=1;}

  /* fill d_ with diagonal values and take reciprocal */
  for(int j=0;j<ny;j++){
    for(int i=0;i<nx;i++){
      idx = j*nx+i;
      d_[idx] = 1;
      d_[idx] = 1/d_[idx];
    }
  }

  ierr = VecCreateSeq(PETSC_COMM_SELF, N, &b); CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, N, &u); CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, N, &x); CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, N, &d); CHKERRQ(ierr);
  ierr = VecPlaceArray(u,u_); CHKERRQ(ierr);
  ierr = VecPlaceArray(x,x_); CHKERRQ(ierr);
  ierr = VecPlaceArray(d,d_); CHKERRQ(ierr);

  /*------------------ MATRICES -----------------------*/
  ierr = MatCreateShell(PETSC_COMM_SELF,N,N,N,N,(void*)&ctx,&A); CHKERRQ(ierr);
  ierr = MatShellSetOperation(A,MATOP_MULT,(void(*)(void)) multA);
  ierr = MatShellSetOperation(A,MATOP_MULT_TRANSPOSE,(void(*)(void)) multA);

  ierr = MatCreateShell(PETSC_COMM_SELF,N,N,N,N,(void*)&ctx,&F); CHKERRQ(ierr);
  ierr = MatShellSetOperation(F,MATOP_MULT,(void(*)(void)) multF);

  /*--------------FFTW Plan--------------*/

  printf("u\n");
  for(int i=0;i<N;i++) printf("%f\t",u_[i]); printf("\n");
  ierr = MatMult(F,u,x); CHKERRQ(ierr);
  printf("x\n");
  for(int i=0;i<N;i++) printf("%f\t",x_[i]); printf("\n");

  ierr = VecDestroy(&u);CHKERRQ(ierr); ierr = VecDestroy(&b);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr); ierr = MatDestroy(&F);CHKERRQ(ierr);
  free(u_);
  free(x_);
  ierr = PetscFinalize();
  return ierr;
}
