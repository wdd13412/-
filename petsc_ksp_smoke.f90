program petsc_ksp_smoke
#include "petsc/finclude/petscksp.h"
  use petscksp
  implicit none

  PetscErrorCode :: ierr
  Mat            :: A
  Vec            :: x, b
  KSP            :: ksp
  PC             :: pc
  PetscInt       :: n, i, Ii, one
  PetscScalar    :: v_diag, v_off
  PetscReal      :: rnorm

  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
  n = 8
  one = 1
  v_diag = 2.0
  v_off  = -1.0

  call MatCreateSeqAIJ(PETSC_COMM_SELF, n, n, 3, PETSC_NULL_INTEGER, A, ierr)
  do i = 0, n-1
    Ii = i
    if (i > 0) call MatSetValue(A, Ii, Ii-1, v_off, INSERT_VALUES, ierr)
    call MatSetValue(A, Ii, Ii, v_diag, INSERT_VALUES, ierr)
    if (i < n-1) call MatSetValue(A, Ii, Ii+1, v_off, INSERT_VALUES, ierr)
  end do
  call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr)
  call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr)

  call VecCreateSeq(PETSC_COMM_SELF, n, b, ierr)
  call VecDuplicate(b, x, ierr)
  call VecSet(b, 1.0d0, ierr)
  call VecSet(x, 0.0d0, ierr)

  call KSPCreate(PETSC_COMM_SELF, ksp, ierr)
  call KSPSetOperators(ksp, A, A, ierr)
  call KSPSetType(ksp, KSPGMRES, ierr)
  call KSPGetPC(ksp, pc, ierr)
  call PCSetType(pc, PCILU, ierr)
  call KSPSetTolerances(ksp, 1.0d-12, PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL, 200, ierr)
  call KSPSetFromOptions(ksp, ierr)
  call KSPSolve(ksp, b, x, ierr)
  call KSPGetResidualNorm(ksp, rnorm, ierr)

  write(*,*) '[PETSc-smoke] residual norm = ', rnorm

  call KSPDestroy(ksp, ierr)
  call VecDestroy(x, ierr)
  call VecDestroy(b, ierr)
  call MatDestroy(A, ierr)
  call PetscFinalize(ierr)
end program petsc_ksp_smoke
