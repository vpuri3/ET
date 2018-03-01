static char help[] = "Solves a tridiagonal linear system with KSP.\n\n";

#include <petscksp.h>
#include <math.h>

typedef struct {
  PetscInt  nx, ny;
} Ctx;

PetscErrorCode multA(Mat B, Vec a, Vec b){
  PetscInt idx;
  Ctx *aa = NULL;
  PetscErrorCode ierr;
  ierr = MatShellGetContext(B,&aa);
  PetscInt nx = aa->nx; PetscInt ny = aa->ny; 
  PetscReal dxinv = (double) nx+1;
  PetscReal dyinv = (double) ny+1;
  PetscReal * a_, * b_;
  ierr = VecGetArray(a,&a_); CHKERRQ(ierr);
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
  VecRestoreArray(a,&a_);
  VecRestoreArray(b,&b_);
  return 0;
}

int main(int argc,char **args)
{
  PetscErrorCode ierr;
  Vec         x, b, u; /* approx solution, RHS, exact solution */
  Mat         A, M, Vx, Vy; /* shell matrix, shell preconditioner matrix, matrix of eigenvectors (column) */
  KSP         ksp;     /* linear solver context */
  PC          pc;      /* preconditioner context */
  PetscReal   norm, xx, yy;
  PetscInt    k, nx=50, ny=50, N=nx*ny, its; /* k is wavenumber */
  PetscBool   nonzeroguess = PETSC_FALSE,changepcside = PETSC_FALSE;
  Ctx         actx;

  ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,NULL,"-nx",&nx,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-ny",&ny,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-nonzero_guess",&nonzeroguess,NULL);CHKERRQ(ierr);

  /*===== QUESTIONS ================
  1. Sequential AIJ sparse matrices for matrix of eigen vectors of matrix A.
  2. vector for matrix of eigenvalues of matrix A.
  5. look into MatSetOption(), KSPSetOption(), VecSetOption().
  -Does Petsc support boundary nodes/ some inbuilt way to add in boundary conditions during matrix multiplication?
  -Is it possible to have different Boundary Conditions on different edges?
  */

  /* Giving problem parameters to context */
  actx.nx = nx; actx.ny = ny;

  /*------------------ MATRIX -----------------------
     When using MatCreate(), the matrix format can
     be specified at runtime.
     Performance tuning note:  For problems of substantial size,
     preallocation of matrix memory is crucial for attaining good
     performance. See the matrix chapter of the users manual for details.
  */

  ierr = MatCreate(PETSC_COMM_SELF,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);

  /* Assemble Matrix V*/
  /*
  value[0] = -1.0; value[1] = 2.0; value[2] = -1.0;
  for (int i=1; i<n-1; i++) {
    col[0] = i-1; col[1] = i; col[2] = i+1;
    ierr   = MatSetValues(V,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);
  }
  i    = n - 1; col[0] = n - 2; col[1] = n - 1;
  ierr = MatSetValues(V,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
  i    = 0; col[0] = 0; col[1] = 1; value[0] = 2.0; value[1] = -1.0;
  ierr = MatSetValues(V,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(V,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(V,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  */

  /* make shell matrix A*/
  ierr = MatCreateShell(PETSC_COMM_SELF,N,N,N,N,(void*)&actx,&A); CHKERRQ(ierr);
  ierr = MatShellSetOperation(A,MATOP_MULT,(void(*)(void)) multA);
  ierr = MatShellSetOperation(A,MATOP_MULT_TRANSPOSE,(void(*)(void)) multA);
  

  /*-----------VECTORS--------------*/
  ierr = VecCreateSeq(PETSC_COMM_SELF, n, &b); CHKERRQ(ierr);
  ierr = VecSetFromOptions(b);CHKERRQ(ierr);
  ierr = VecDuplicate(b,&x);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) x, "Solution");CHKERRQ(ierr);

  /* array pp --> vector u, array q --> vector x. all malloc allocations are freed at the end */
  double * pp = malloc(n*sizeof(double));
  double * q = malloc(n*sizeof(double));
  double nn = n;  double dx = 1/(nn+1);
  for(int i=0;i<n;i++){
    double xx = (i+1)*dx;
    q[i] = xx*xx;
    pp[i] = 0;
  }

  ierr = VecPlaceArray(x,pp);
  ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1,n,q,&u); CHKERRQ(ierr);
  ierr = MatMult(B,u,b);CHKERRQ(ierr);
  //ierr = multiply(B,u,b);

  /* Create KSP object */
  ierr = KSPCreate(PETSC_COMM_SELF,&ksp);CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,B,A);CHKERRQ(ierr);

  /* Test if you can change the KSP side and type after they have been previously set */
  ierr = PetscOptionsGetBool(NULL,NULL,"-change_pc_side",&changepcside,NULL);CHKERRQ(ierr);
  if (changepcside) {
    ierr = KSPSetUp(ksp);CHKERRQ(ierr);
    ierr = PetscOptionsInsertString(NULL,"-ksp_norm_type unpreconditioned -ksp_pc_side right");CHKERRQ(ierr);
    ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
    ierr = KSPSetUp(ksp);CHKERRQ(ierr);
  }

  /*
     Set linear solver defaults for this problem (optional).
     - By extracting the KSP and PC contexts from the KSP context,
       we can then directly call any KSP and PC routines to set
       various options.
     - The following four statements are optional; all of these
       parameters could alternatively be specified at runtime via
       KSPSetFromOptions();
  */
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  ierr = PCSetType(pc,PCJACOBI);CHKERRQ(ierr);
  ierr = KSPSetTolerances(ksp,1.e-5,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);

  /*
    Set runtime options, e.g.,
        -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
    These options will override those specified above as long as
    KSPSetFromOptions() is called _after_ any other customization
    routines.
  */
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);

  if (nonzeroguess) {
    PetscScalar p = .5;
    ierr = VecSet(x,p);CHKERRQ(ierr);
    ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);CHKERRQ(ierr);
  }

  /* Solve linear system */
  ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);

  /* View solver info (can instead use -ksp_view) */
  ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);


  /* 
  printf("x\t");
  ierr = VecView(x, PETSC_VIEWER_STDOUT_SELF);
  */

  /* Check the error*/
  ierr = VecAXPY(x,-1.0,u);CHKERRQ(ierr);
  ierr = VecNorm(x,NORM_2,&norm);CHKERRQ(ierr);
  ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF,"Norm of error %g, Iterations %D\n",(double)norm,its);CHKERRQ(ierr);

  /*
     Free work space.  All PETSc objects should be destroyed when they
     are no longer needed.
  */
  ierr = VecDestroy(&x);CHKERRQ(ierr); ierr = VecDestroy(&u);CHKERRQ(ierr);
  ierr = VecDestroy(&b);CHKERRQ(ierr); ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);

  free(pp);
  free(q);
  /*
     Always call PetscFinalize() before exiting a program.  This routine
       - finalizes the PETSc libraries as well as MPI
       - provides summary and diagnostic information if certain runtime
         options are chosen (e.g., -log_view).
  */
  ierr = PetscFinalize();
  return ierr;
}
