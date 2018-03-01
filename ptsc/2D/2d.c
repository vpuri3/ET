static char help[] = "Solves a tridiagonal linear system with KSP.\n\n";

#include <petscksp.h>
#include <math.h>

typedef struct {
  PetscInt  nx, ny;
} Ctx;

PetscErrorCode multA(Mat A, Vec a, Vec b){
  PetscInt idx;
  Ctx *aa = NULL;
  PetscErrorCode ierr;
  ierr = MatShellGetContext(A,&aa);
  PetscInt nx = aa->nx; PetscInt ny = aa->ny; 
  PetscReal dxinv = (double) nx+1;
  PetscReal dyinv = (double) ny+1;
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

int main(int argc,char **args)
{
  PetscErrorCode ierr;
  Vec         u, b; 
  Mat         A; /* Second order Laplace Operator */
  PetscInt    nx=10, ny=10, N=nx*ny; /* k is wavenumber */
  Ctx         actx;

  ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,NULL,"-nx",&nx,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-ny",&ny,NULL);CHKERRQ(ierr);

  /* putting in context information */
  actx.nx = nx; actx.ny = ny;

  /*------------------ MATRIX -----------------------
    below command may not be necessary. Not sure rn.
  ierr = MatCreate(PETSC_COMM_SELF,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);
  */

  ierr = MatCreateShell(PETSC_COMM_SELF,N,N,N,N,(void*)&actx,&A); CHKERRQ(ierr);
  ierr = MatShellSetOperation(A,MATOP_MULT,(void(*)(void)) multA);
  ierr = MatShellSetOperation(A,MATOP_MULT_TRANSPOSE,(void(*)(void)) multA);

  /*----CREATING VECTORS--------------*/
  PetscReal * u_ = malloc(N*sizeof(PetscReal));
  //PetscReal * b_ = malloc(N*sizeof(PetscReal));
  for(int i=0;i<N;i++) u_[i] = 1; //b_[i] = 0;

  ierr = VecCreateSeq(PETSC_COMM_SELF, N, &b); CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, N, &u); CHKERRQ(ierr);
  //ierr = VecDuplicate(u,&b);CHKERRQ(ierr);
  ierr = VecPlaceArray(u,u_); CHKERRQ(ierr);

  /*-------Doing stuff with vectors-----*/
  ierr = MatMult(A,u,b);CHKERRQ(ierr);

  printf("b\t");
  ierr = VecView(b, PETSC_VIEWER_STDOUT_SELF);

  ierr = VecDestroy(&u);CHKERRQ(ierr); ierr = VecDestroy(&b);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  free(u_);
  ierr = PetscFinalize();
  return ierr;
}
