
static char help[] = "Solves scaled Laplace eqn with variable coefficients in 2D.\n\n";

/* solves cxx*u_xx + cyy*u_yy + c1*u = r. Domain: [0,1]^2. grid size: nx*ny.
   nx = total # of points in x.
   # of interior points in x = nx-2
   dx = 1/(nx-1)

   Turn inhomogeneous problem into a homogeneous problem by subtracting in
   the boundary values first and then zeroing them.
   Key is matrix A preserves boundary elements.
*/

#include <petscksp.h>
#include <math.h>
#include <fftw3.h>
#include "3d.h"

int main(int argc,char **args)
{
  PetscErrorCode ierr;
  Vec            x,u,f,d,cxx,cyy,czz,cave,c1,c1u;  /*sol, true, rhs, diagonals, coeffs*/
  Mat            A, F;                        /*system matrix Operator, preconditioner matrix*/
  PetscInt       idx, nx=50, ny=50, nz=50, N=nx*ny*nz;
  Ctx            ctx;                         /*data structure to pass information around*/
  KSP            ksp;
  PC             pc;
  PetscReal      dx, dy, dz, dxinv, dyinv, dzinv, fft_factor; /* to pass to ctx */

  dxinv = (double) nx-1; dyinv = (double) ny-1; dzinv = (double) nz-1;
  dx = 1/dxinv; dy = 1/dyinv; dz = 1/dzinv; dz = 1/dzinv;
  fft_factor = (double) 0.125/(nz+1)/(ny+1)/(nx+1);

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
  PetscReal * czz_ = malloc(N*sizeof(PetscReal));
  PetscReal * cave_ = malloc(N*sizeof(PetscReal));
  PetscReal * c1_ = malloc(N*sizeof(PetscReal));
  PetscReal * c1u_ = malloc(N*sizeof(PetscReal));

  /*********** Putting some values in **********/
  PetscReal xx=0,yy=0,zz=0, ii=0,jj=0,kk=0;
  for(int k=0;k<nz;k++){
    for(int j=0;j<ny;j++){
      for(int i=0;i<nx;i++){
	idx = k*nx*ny+j*nx+i;
	xx = i*dx; yy = j*dy; ii = i+1; jj = j+1;
	/* eigen values and coefficients */
	d_[idx]  = -2*dxinv*dxinv*(1-cos(ii*M_PI*dx));
	d_[idx] += -2*dyinv*dyinv*(1-cos(jj*M_PI*dy));
	d_[idx] += -2*dzinv*dzinv*(1-cos(kk*M_PI*dz)); /* x eig val + y eig val */
	cxx_[idx] = c_xx(xx,yy,zz);
	cyy_[idx] = c_yy(xx,yy,zz);
	czz_[idx] = c_zz(xx,yy,zz);
	c1_[idx] = c_1(xx,yy,zz);
	/* exact sol, initial guess, rhs */
	u_[idx] = exact(xx,yy,zz);
	x_[idx] = 0;
	f_[idx] = rhs(xx,yy,zz);
	if(Bflag(i,j,k,nx,ny,nz)){
	  x_[idx] = exact(xx,yy,zz); // contains boundary values.
	  f_[idx] = exact(xx,yy,zz);
	}
      }
    }
  }

  /*----------------Putting things in context-------------*/
  ctx.nx = nx; ctx.ny = ny; ctx.nz = nz; ctx.N = N;
  ctx.dx = dx; ctx.dy = dy; ctx.dz = dz;
  ctx.dxinv = dxinv; ctx.dyinv = dyinv; ctx.dzinv = dzinv;
  ctx.fft_factor = fft_factor;
  ctx.dptr = &d; ctx.cxxptr = &cxx; ctx.cyyptr = &cyy; ctx.caveptr = &cave;
  ctx.c1ptr = &c1; ctx.c1uptr = &c1u;
  ctx.Fptr = &F;
  ctx.cxx_ = cxx_; ctx.cyy_ = cyy_; ctx.czz_ = czz_;
  ctx.c1_ = c1_; ctx.c1u_ = c1u_;

  /*----------------CREATING VECTORS--------------------*/
  ierr = VecCreateSeq(PETSC_COMM_SELF, N, &x); CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, N, &u); CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, N, &f); CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, N, &d); CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, N, &cxx); CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, N, &cyy); CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, N, &czz); CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, N, &cave); CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, N, &c1); CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, N, &c1u); CHKERRQ(ierr);
  ierr = VecPlaceArray(x,x_); CHKERRQ(ierr);
  ierr = VecPlaceArray(u,u_); CHKERRQ(ierr);
  ierr = VecPlaceArray(f,f_); CHKERRQ(ierr);
  ierr = VecPlaceArray(d,d_); CHKERRQ(ierr);
  ierr = VecPlaceArray(cxx,cxx_); CHKERRQ(ierr);
  ierr = VecPlaceArray(cyy,cyy_); CHKERRQ(ierr);
  ierr = VecPlaceArray(czz,czz_); CHKERRQ(ierr);
  ierr = VecPlaceArray(cave,cave_); CHKERRQ(ierr);
  ierr = VecPlaceArray(c1,c1_); CHKERRQ(ierr);
  ierr = VecPlaceArray(c1u,c1u_); CHKERRQ(ierr);

  ierr = VecReciprocal(d);
  ierr = VecCopy(cxx,cave);
  ierr = VecAXPBY(cave,1/3,1/3,cyy);
  ierr = VecAXPBY(cave,1,1/3,czz);
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
  if(precon){
    ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
    ierr = PCSetType(pc,PCSHELL); CHKERRQ(ierr);
    
    ierr = PCShellSetApply(pc,precondition);
    ierr = PCShellSetContext(pc,&ctx);
  }
  ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);
  ierr = KSPSetType(ksp,KSPBCGS);CHKERRQ(ierr);
  ierr = KSPSetTolerances(ksp,1e-14,1e-14,PETSC_DEFAULT,PETSC_DEFAULT);
  ierr = KSPSetNormType(ksp,KSP_NORM_UNPRECONDITIONED);
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);

  /*------------------Solving-----------------------*/
  Vec v1, v2; /* general purpose vector */
  ierr = VecCreateSeq(PETSC_COMM_SELF, N, &v1); CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, N, &v2); CHKERRQ(ierr);

  ierr = MatMult(A,x,v1);
  ierr = VecAXPY(f,-1,v1); /* subtracting boundary values. f = f - Ax */

  ierr = VecCopy(x,v1);
  ierr = VecZeroEntries(x);
  ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);
  ierr = KSPSolve(ksp,f,x);
  ierr = VecAXPY(x,1,v1); /* add in boundary values */
  ierr = MatMult(A,v1,v2);
  ierr = VecAXPY(f,1,v2);

  /* check if |Au - f| */
  double e;
  ierr = MatMult(A,u,v1);
  ierr = VecAXPY(v1,-1,f);
  ierr = VecNorm(v1,NORM_2,&e);
  ierr = PetscPrintf(PETSC_COMM_SELF,"2 Norm of |Au-f| is %g.\n",(double)e);CHKERRQ(ierr);
  /* get residual */
  ierr = VecZeroEntries(v1);
  ierr = MatMult(A,x,v1);
  ierr = VecAXPY(v1,-1,f);
  ierr = VecNorm(v1,NORM_2,&e);
  ierr = PetscPrintf(PETSC_COMM_SELF,"2 Norm of |Ax-f| is %g.\n",(double)e);CHKERRQ(ierr);
  /* get error */
  ierr = VecAXPY(x,-1,u);
  ierr = VecNorm(x,NORM_INFINITY,&e);
  ierr = PetscPrintf(PETSC_COMM_SELF,"Inf Norm of |u-x| is %g.\n",(double)e);CHKERRQ(ierr);

  /*-------------------End Credits---------------------*/
  ierr = VecDestroy(&x);CHKERRQ(ierr); ierr = VecDestroy(&u);CHKERRQ(ierr);
  ierr = VecDestroy(&f);CHKERRQ(ierr); ierr = VecDestroy(&d);CHKERRQ(ierr);
  //ierr = VecDestry(&v1);CHKERRQ(ierr); ierr = VecDestry(&v1);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr); ierr = MatDestroy(&F);CHKERRQ(ierr);
  ierr = KSPDestroy(&ksp);
  free(x_); free(u_);
  free(f_); free(d_);
  ierr = PetscFinalize();
  return ierr;
}
