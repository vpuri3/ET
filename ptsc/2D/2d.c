static char help[] = "Solves a tridiagonal linear system with KSP.\n\n";

#include <petscksp.h>
#include <math.h>

typedef struct {
  PetscInt  nx, ny;
  PetscReal dx, dy, dxinv, dyinv;
} Ctx;

PetscErrorCode multA(Mat A, Vec a, Vec b){
  PetscInt idx;
  Ctx *aa = NULL;
  PetscErrorCode ierr;
  ierr = MatShellGetContext(A,&aa);
  PetscInt nx = aa->nx; PetscInt ny = aa->ny; 
  double dxinv = aa->dxinv;
  double dyinv = aa->dyinv;
  const PetscReal * a_ = NULL; PetscReal * b_ = NULL;
  ierr = VecGetArrayRead(a,&a_); CHKERRQ(ierr);
  ierr = VecGetArray(b,&b_); CHKERRQ(ierr);
  PetscReal xm1, xp1, ym1, yp1;
  for(int j=0;j<ny;j++){
    for(int i=0;i<nx;i++){
      idx = j*ny + i;
      if (i == 0) xm1 = 0;
      else xm1 = a_[idx-1];
      if(i == nx-1) xp1 = 0;
      else xp1 = a_[idx+1];
      if (j == 0) ym1 = 0;
      else ym1 = a_[idx-nx];
      if(j == ny-1) yp1 = 0;
      else yp1 = a_[idx+nx];
      b_[idx] = dxinv*dxinv*(-xm1+2*a_[idx]-xp1) + dyinv*dyinv*(-ym1+2*a_[idx]-yp1);
    }
  }
  VecRestoreArrayRead(a,&a_);
  VecRestoreArray(b,&b_);
  return 0;
}

/* Learn to use FFTW with PETSc fore preconditioner matrix. */

int main(int argc,char **args)
{
  PetscErrorCode ierr;
  Vec         u, b; 
  Mat         A, Vx, Vy; /* Second order Laplace Operator */
  PetscInt    nx=10, ny=10, N=nx*ny;
  Ctx         ctx;

  ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,NULL,"-nx",&nx,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-ny",&ny,NULL);CHKERRQ(ierr);

  /*----------computing constants and stuff----------*/
  ctx.nx = nx; ctx.ny = ny;
  ctx.dxinv = (double) nx+1; ctx.dyinv = (double) ny+1;
  ctx.dx = 1/ctx.dxinv; ctx.dy = 1/ctx.dyinv;

  /*------------------ MATRICES -----------------------*/
  ierr = MatCreateShell(PETSC_COMM_SELF,N,N,N,N,(void*)&ctx,&A); CHKERRQ(ierr);
  ierr = MatShellSetOperation(A,MATOP_MULT,(void(*)(void)) multA);
  ierr = MatShellSetOperation(A,MATOP_MULT_TRANSPOSE,(void(*)(void)) multA);

  ierr = MatCreateSeqDense(PETSC_COMM_SELF,nx,nx,NULL,&Vx); CHKERRQ(ierr);

  /*----------------CREATING VECTORS--------------*/
  PetscReal * u_ = malloc(N*sizeof(PetscReal));
  for(int i=0;i<N;i++) u_[i] = 1;

  ierr = VecCreateSeq(PETSC_COMM_SELF, N, &b); CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, N, &u); CHKERRQ(ierr);
  ierr = VecPlaceArray(u,u_); CHKERRQ(ierr);

  ierr = VecDestroy(&u);CHKERRQ(ierr); ierr = VecDestroy(&b);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  free(u_);
  ierr = PetscFinalize();
  return ierr;
}
