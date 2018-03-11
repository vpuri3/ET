
static char help[] = "Solves scaled Laplace eqn with variable coefficients in 2D.\n\n";

/* solves cxx*u_xx + cyy*u_yy + c1*u = r. Domain: [0,1]^2. grid size: nx*ny */


#include <petscksp.h>
#include <math.h>
#include <fftw3.h>
#include "2d.h"

int main(int argc,char **args)
{
  PetscErrorCode ierr;
  Vec            x, u, f, d, cxx, cyy, cave, c1, temp;  /*unknown, true, rhs, diagonals, coeffs*/
  Mat            A, F;                        /*second order Laplace Operator, fft matrix*/
  PetscInt       idx, nx=500, ny=500, N=nx*ny;
  Ctx            ctx;                         /*data structure to pass information around*/
  KSP            ksp;
  PC             pc;
  PetscReal      dx, dy, dxinv, dyinv, fft_factor; /* to pass to ctx */

  dxinv = (double) nx+1; dyinv = (double) ny+1;
  dx = 1/dxinv; dy = 1/dyinv;
  fft_factor = (double) 0.25/(ny+1)/(nx+1);

  ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,NULL,"-nx",&nx,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-ny",&ny,NULL);CHKERRQ(ierr);

  /*----------------Allocating vectors-------------*/
  PetscReal * x_ = malloc(N*sizeof(PetscReal));
  PetscReal * u_ = malloc(N*sizeof(PetscReal));
  PetscReal * f_ = malloc(N*sizeof(PetscReal));
  PetscReal * d_ = malloc(N*sizeof(PetscReal));
  PetscReal * cxx_ = malloc(N*sizeof(PetscReal));
  PetscReal * cyy_ = malloc(N*sizeof(PetscReal));
  PetscReal * cave_ = malloc(N*sizeof(PetscReal));
  PetscReal * c1_ = malloc(N*sizeof(PetscReal));
  PetscReal * temp_ = malloc(N*sizeof(PetscReal));

  PetscReal xx=0,yy=0, ii=0, jj=0;
  for(int j=0;j<ny;j++){ for(int i=0;i<nx;i++){
      idx = j*nx+i;
      xx = i*dx; yy = j*dy; ii = i+1; jj = j+1;
      u_[idx] = xx*xx + exp(yy);
      x_[idx] = 0.1;
      f_[idx] = 0;
      d_[idx] = 2*dxinv*dxinv*(1-cos(ii*M_PI*dx));
      d_[idx] += 2*dyinv*dyinv*(1-cos(jj*M_PI*dy)); /* x eig val + y eig val */
      cxx_[idx] = exp(xx+yy);
      cyy_[idx] = yy*yy*yy+1;
      c1_[idx] = xx + exp(yy*yy);
    }
  }

  /*----------------Putting things in context-------------*/
  ctx.nx = nx; ctx.ny = ny; ctx.N = N;
  ctx.dx = dx; ctx.dy = dy; ctx.dxinv = dxinv; ctx.dyinv = dyinv;
  ctx.fft_factor = fft_factor;
  ctx.dptr = &d; ctx.cxxptr = &cxx; ctx.cyyptr = &cyy; ctx.caveptr = &cave;
  ctx.c1ptr = &c1; ctx.tempptr = &temp;
  ctx.Fptr = &F;
  ctx.cxx_ = cxx_; ctx.cyy_ = cyy_;

  /*----------------CREATING VECTORS--------------------*/
  ierr = VecCreateSeq(PETSC_COMM_SELF, N, &x); CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, N, &u); CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, N, &f); CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, N, &d); CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, N, &cxx); CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, N, &cyy); CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, N, &cave); CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, N, &c1); CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, N, &temp); CHKERRQ(ierr);
  ierr = VecPlaceArray(x,x_); CHKERRQ(ierr);
  ierr = VecPlaceArray(u,u_); CHKERRQ(ierr);
  ierr = VecPlaceArray(f,f_); CHKERRQ(ierr);
  ierr = VecPlaceArray(d,d_); CHKERRQ(ierr);
  ierr = VecPlaceArray(cxx,cxx_); CHKERRQ(ierr);
  ierr = VecPlaceArray(cyy,cyy_); CHKERRQ(ierr);
  ierr = VecPlaceArray(cave,cave_); CHKERRQ(ierr);
  ierr = VecPlaceArray(c1,c1_); CHKERRQ(ierr);
  ierr = VecPlaceArray(temp,temp_); CHKERRQ(ierr);

  ierr = VecReciprocal(d);
  ierr = VecCopy(cxx,cave);
  ierr = VecAXPBY(cave,0.5,0.5,cyy);
  ierr = VecReciprocal(cave);
  /*------------------ MATRICES -----------------------*/
  ierr = MatCreateShell(PETSC_COMM_SELF,N,N,N,N,(void*)&ctx,&A); CHKERRQ(ierr);
  ierr = MatShellSetOperation(A,MATOP_MULT,(void(*)(void)) multA);
  ierr = MatShellSetOperation(A,MATOP_MULT_TRANSPOSE,(void(*)(void)) multATransp);

  ierr = MatCreateShell(PETSC_COMM_SELF,N,N,N,N,(void*)&ctx,&F); CHKERRQ(ierr);
  ierr = MatShellSetOperation(F,MATOP_MULT,(void(*)(void)) multF);
  ierr = MatShellSetOperation(F,MATOP_MULT_TRANSPOSE,(void(*)(void)) multFTransp);

  /*----------------------KSP-------------------------*/
  ierr = KSPCreate(PETSC_COMM_SELF,&ksp);CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);

  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  ierr = PCSetType(pc,PCSHELL); CHKERRQ(ierr); /* look up PCShellSetApply() when using PCSHELL */

  ierr = PCShellSetApply(pc,precondition);
  ierr = PCShellSetContext(pc,&ctx);

  ierr = KSPSetType(ksp,KSPBCGS);CHKERRQ(ierr);
  ierr = KSPSetTolerances(ksp,1e-12,1e-12,PETSC_DEFAULT,PETSC_DEFAULT);
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);

  /*------------------Solving-----------------------*/
  ierr = MatMult(A,u,f);
  ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);
  ierr = KSPSolve(ksp,f,x);
  /* get error */
  ierr = VecAXPY(x,-1,u);
  double e;
  ierr = VecNorm(x,NORM_2,&e);
  ierr = PetscPrintf(PETSC_COMM_SELF,"Norm of error is %g.",(double)e);CHKERRQ(ierr);

  /*-------------------End Credits---------------------*/
  ierr = VecDestroy(&x);CHKERRQ(ierr); ierr = VecDestroy(&u);CHKERRQ(ierr);
  ierr = VecDestroy(&f);CHKERRQ(ierr); ierr = VecDestroy(&d);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr); ierr = MatDestroy(&F);CHKERRQ(ierr);
  ierr = KSPDestroy(&ksp);
  free(x_); free(u_);
  free(f_); free(d_);
  ierr = PetscFinalize();
  return ierr;
}
