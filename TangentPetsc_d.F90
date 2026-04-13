MODULE TANGENT_PETSC_MODULE
#include <petsc/finclude/petscksp.h>
  USE petscksp
  USE BUFLOWMODULE_DIFF, ONLY: TANGENT_MATVEC
  IMPLICIT NONE

  LOGICAL, SAVE :: petsc_initialized_local = .FALSE.
  INTEGER(kind=8), SAVE :: tangent_n_ctx = 0_8
  INTEGER(kind=8), SAVE :: tangent_matmult_calls = 0_8
  INTEGER(kind=8), SAVE :: tangent_idx_n = 0_8
  REAL(kind=8), ALLOCATABLE, SAVE :: tangent_data_ctx(:, :)
  REAL(kind=8), ALLOCATABLE, SAVE :: tangent_cell_ctx(:, :)
  PetscInt, ALLOCATABLE, SAVE :: tangent_idx(:)

CONTAINS

  SUBROUTINE TANGENT_MATVEC_WRAP(nloc, xin, yout)
    IMPLICIT NONE
    INTEGER(kind=8), INTENT(IN) :: nloc
    REAL(kind=8), INTENT(IN) :: xin(nloc)
    REAL(kind=8), INTENT(OUT) :: yout(nloc)
    CALL TANGENT_MATVEC(tangent_data_ctx, tangent_cell_ctx, nloc, xin, yout)
  END SUBROUTINE TANGENT_MATVEC_WRAP

  SUBROUTINE TANGENT_PETSC_SET_CONTEXT(data_4d137, cellprimitives, n)
    IMPLICIT NONE
    REAL(kind=8), INTENT(IN) :: data_4d137(1, 4)
    REAL(kind=8), INTENT(IN) :: cellprimitives(:, :)
    INTEGER(kind=8), INTENT(IN) :: n
    INTEGER(kind=8) :: nc, i

    nc = SIZE(cellprimitives, 1)
    tangent_n_ctx = n

    IF (ALLOCATED(tangent_data_ctx)) DEALLOCATE(tangent_data_ctx)
    IF (ALLOCATED(tangent_cell_ctx)) DEALLOCATE(tangent_cell_ctx)
    ALLOCATE(tangent_data_ctx(1, 4), tangent_cell_ctx(nc, 5))
    tangent_data_ctx = data_4d137
    tangent_cell_ctx = cellprimitives(:, 1:5)
    IF (ALLOCATED(tangent_idx)) DEALLOCATE(tangent_idx)
    ALLOCATE(tangent_idx(n))
    DO i = 1_8, n
      tangent_idx(i) = i - 1_8
    END DO
    tangent_idx_n = n
  END SUBROUTINE TANGENT_PETSC_SET_CONTEXT

  SUBROUTINE TANGENT_PETSC_MATMULT_IMPL(A, X, Y, ierr)
    IMPLICIT NONE
    Mat :: A
    Vec :: X, Y
    PetscErrorCode :: ierr
    PetscScalar, ALLOCATABLE :: xvals(:), yvals(:)
    REAL(kind=8), ALLOCATABLE :: xin(:), yout(:)
    INTEGER(kind=8) :: nloc, i
    LOGICAL :: had_bad
    REAL(kind=8) :: xin_norm, yout_norm

    ierr = 0
    nloc = tangent_n_ctx
    IF (nloc <= 0_8) RETURN

    IF (tangent_idx_n /= nloc .OR. .NOT. ALLOCATED(tangent_idx)) THEN
      ierr = 1
      RETURN
    END IF
    ALLOCATE(xin(nloc), yout(nloc), xvals(nloc), yvals(nloc))

    CALL VecGetValues(X, INT(nloc, kind=KIND(1)), tangent_idx, xvals, ierr)
    IF (ierr /= 0) THEN
      DEALLOCATE(xin, yout, xvals, yvals)
      RETURN
    END IF
    xin = xvals

    had_bad = .FALSE.
    DO i = 1_8, nloc
      IF (xin(i) /= xin(i) .OR. ABS(xin(i)) > 1.0d300) THEN
        xin(i) = 0.0_8
        had_bad = .TRUE.
      END IF
    END DO

    CALL TANGENT_MATVEC_WRAP(nloc, xin, yout)
    IF (tangent_matmult_calls == 0_8) THEN
      PRINT *, '[PETSC-MATMULT-DBG] raw yout(1:5)=', yout(1:MIN(5_8, nloc))
    END IF

    DO i = 1_8, nloc
      IF (yout(i) /= yout(i) .OR. ABS(yout(i)) > 1.0d300) THEN
        IF (.NOT. had_bad .AND. tangent_matmult_calls <= 2_8) THEN
          PRINT *, '[PETSC-MATMULT-DBG] first bad yout idx=', i, ' value=', yout(i)
        END IF
        yout(i) = 0.0_8
        had_bad = .TRUE.
      END IF
    END DO

    yvals = yout
    CALL VecSetValues(Y, INT(nloc, kind=KIND(1)), tangent_idx, yvals, INSERT_VALUES, ierr)
    CALL VecAssemblyBegin(Y, ierr)
    CALL VecAssemblyEnd(Y, ierr)

    tangent_matmult_calls = tangent_matmult_calls + 1_8
    IF (tangent_matmult_calls <= 3_8) THEN
      xin_norm = SQRT(SUM(xin*xin))
      yout_norm = SQRT(SUM(yout*yout))
      PRINT *, '[PETSC-MATMULT] call=', tangent_matmult_calls, ' ||x||=', xin_norm, ' ||Ax||=', yout_norm, ' sanitized=', had_bad
    END IF

    DEALLOCATE(xin, yout, xvals, yvals)
  END SUBROUTINE TANGENT_PETSC_MATMULT_IMPL

  SUBROUTINE SOLVE_TANGENT_PETSC_GMRES(data_4d137, cellprimitives, n, b, x, tol, maxiter, m_restart, info)
    IMPLICIT NONE
    REAL(kind=8), INTENT(IN) :: data_4d137(1, 4), cellprimitives(:, :), b(n), tol
    REAL(kind=8), INTENT(INOUT) :: x(n)
    INTEGER(kind=8), INTENT(IN) :: n, maxiter, m_restart
    INTEGER, INTENT(OUT) :: info

    Mat :: A
    Vec :: vb, vx, vtmp
    KSP :: ksp
    PC :: pc
    KSPConvergedReason :: reason
    PetscErrorCode :: ierr
    PetscInt :: n_petsc, max_it_petsc, restart_petsc, its, i
    PetscReal :: rnorm, bnorm, xnorm0, axnorm0, axnormb
    PetscInt, ALLOCATABLE :: idx(:)
    PetscScalar, ALLOCATABLE :: vals(:)
    EXTERNAL TANGENT_PETSC_MATMULT

    info = 2
    IF (n <= 0_8) RETURN

    CALL TANGENT_PETSC_SET_CONTEXT(data_4d137, cellprimitives, n)

    IF (.NOT. petsc_initialized_local) THEN
      CALL PetscInitialize(PETSC_NULL_CHARACTER, ierr)
      IF (ierr /= 0) THEN
        PRINT *, '[PETSC] initialize failed, ierr=', ierr
        RETURN
      END IF
      petsc_initialized_local = .TRUE.
    END IF

    n_petsc = INT(n, KIND=KIND(n_petsc))
    max_it_petsc = INT(MAX(1_8, maxiter) * MAX(1_8, m_restart), KIND=KIND(max_it_petsc))
    restart_petsc = INT(MAX(5_8, m_restart), KIND=KIND(restart_petsc))

    CALL MatCreateShell(PETSC_COMM_SELF, n_petsc, n_petsc, n_petsc, n_petsc, PETSC_NULL_INTEGER, A, ierr)
    IF (ierr /= 0) THEN
      PRINT *, '[PETSC] MatCreateShell failed, ierr=', ierr
      RETURN
    END IF
    CALL MatShellSetOperation(A, MATOP_MULT, TANGENT_PETSC_MATMULT, ierr)
    IF (ierr /= 0) THEN
      PRINT *, '[PETSC] MatShellSetOperation failed, ierr=', ierr
      CALL MatDestroy(A, ierr)
      RETURN
    END IF

    CALL VecCreateSeq(PETSC_COMM_SELF, n_petsc, vb, ierr)
    CALL VecDuplicate(vb, vx, ierr)

    ALLOCATE(idx(n_petsc), vals(n_petsc))
    DO i = 1, n_petsc
      idx(i) = i - 1
    END DO
    vals = b(1:n)
    CALL VecSetValues(vb, n_petsc, idx, vals, INSERT_VALUES, ierr)
    IF (ierr /= 0) THEN
      PRINT *, '[PETSC] VecSetValues(vb) failed, ierr=', ierr
      CALL VecDestroy(vb, ierr); CALL VecDestroy(vx, ierr); CALL MatDestroy(A, ierr)
      RETURN
    END IF
    CALL VecAssemblyBegin(vb, ierr)
    CALL VecAssemblyEnd(vb, ierr)
    vals = 0.0d0
    CALL VecSetValues(vx, n_petsc, idx, vals, INSERT_VALUES, ierr)
    IF (ierr /= 0) THEN
      PRINT *, '[PETSC] VecSetValues(vx) failed, ierr=', ierr
      CALL VecDestroy(vb, ierr); CALL VecDestroy(vx, ierr); CALL MatDestroy(A, ierr)
      RETURN
    END IF
    CALL VecAssemblyBegin(vx, ierr)
    CALL VecAssemblyEnd(vx, ierr)
    CALL VecNorm(vb, NORM_2, bnorm, ierr)
    CALL VecNorm(vx, NORM_2, xnorm0, ierr)
    PRINT *, '[PETSC-GMRES] pre-solve ||b||=', bnorm, ' ||x0||=', xnorm0

    CALL VecDuplicate(vx, vtmp, ierr)
    CALL MatMult(A, vx, vtmp, ierr)
    CALL VecNorm(vtmp, NORM_2, axnorm0, ierr)
    CALL MatMult(A, vb, vtmp, ierr)
    CALL VecNorm(vtmp, NORM_2, axnormb, ierr)
    PRINT *, '[PETSC-GMRES] pre-check ||A*x0||=', axnorm0, ' ||A*b||=', axnormb
    CALL VecDestroy(vtmp, ierr)

    CALL KSPCreate(PETSC_COMM_SELF, ksp, ierr)
    CALL KSPSetOperators(ksp, A, A, ierr)
    CALL KSPSetType(ksp, KSPFGMRES, ierr)
    CALL KSPGMRESSetRestart(ksp, restart_petsc, ierr)
    CALL KSPGetPC(ksp, pc, ierr)
    CALL PCSetType(pc, PCNONE, ierr)
    CALL KSPSetTolerances(ksp, tol, PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL, max_it_petsc, ierr)
    CALL KSPSetFromOptions(ksp, ierr)
    CALL KSPSetUp(ksp, ierr)
    CALL KSPSolve(ksp, vb, vx, ierr)
    IF (ierr /= 0) THEN
      PRINT *, '[PETSC] KSPSolve ierr=', ierr
    END IF
    CALL KSPGetConvergedReason(ksp, reason, ierr)
    CALL KSPGetResidualNorm(ksp, rnorm, ierr)
    CALL KSPGetIterationNumber(ksp, its, ierr)

    vals = 0.0_8
    CALL VecGetValues(vx, n_petsc, idx, vals, ierr)
    x(1:n) = vals(1:n_petsc)

    bnorm = SQRT(SUM(b*b))
    PRINT *, '[PETSC-GMRES] reason=', reason, ' its=', its, ' true||r||/||b||=', rnorm / MAX(1.0d-30, bnorm)

    IF (reason > 0) THEN
      info = 0
    ELSE
      info = 1
    END IF

    CALL KSPDestroy(ksp, ierr)
    CALL VecDestroy(vb, ierr)
    CALL VecDestroy(vx, ierr)
    CALL MatDestroy(A, ierr)
    DEALLOCATE(idx, vals)
  END SUBROUTINE SOLVE_TANGENT_PETSC_GMRES

  SUBROUTINE TANGENT_PETSC_FINALIZE()
    IMPLICIT NONE
    PetscErrorCode :: ierr
    IF (petsc_initialized_local) THEN
      CALL PetscFinalize(ierr)
      petsc_initialized_local = .FALSE.
    END IF
    IF (ALLOCATED(tangent_data_ctx)) DEALLOCATE(tangent_data_ctx)
    IF (ALLOCATED(tangent_cell_ctx)) DEALLOCATE(tangent_cell_ctx)
    IF (ALLOCATED(tangent_idx)) DEALLOCATE(tangent_idx)
    tangent_idx_n = 0_8
    tangent_n_ctx = 0_8
  END SUBROUTINE TANGENT_PETSC_FINALIZE

END MODULE TANGENT_PETSC_MODULE

SUBROUTINE TANGENT_PETSC_MATMULT(A, X, Y, ierr)
  USE petscksp
  USE TANGENT_PETSC_MODULE, ONLY: TANGENT_PETSC_MATMULT_IMPL
  IMPLICIT NONE
  Mat :: A
  Vec :: X, Y
  PetscErrorCode :: ierr
  CALL TANGENT_PETSC_MATMULT_IMPL(A, X, Y, ierr)
END SUBROUTINE TANGENT_PETSC_MATMULT
