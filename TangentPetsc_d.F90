MODULE TANGENT_PETSC_MODULE
#include <petsc/finclude/petscksp.h>
  USE petscksp
  USE BUFLOWMODULE_DIFF, ONLY: TANGENT_MATVEC, APPLY_TANGENT_OPERATOR_WORK, &
& APPLY_RESIDUAL_SCALING, pc_row_ptr, pc_col_ind, pc_blk, pc_a0_blk, &
& pc_a0_ready, pc_diag_pos
  USE TYPESMODULE_DIFF, ONLY: point_updated, fluxresiduals_slnd, facefluxes_slnd, &
 & cellfluxes_slnd, cellstate_slnd, cellprimitives_slnd
  IMPLICIT NONE

  LOGICAL, SAVE :: petsc_initialized_local = .FALSE.
  INTEGER(kind=8), SAVE :: tangent_n_ctx = 0_8
  INTEGER(kind=8), SAVE :: tangent_matmult_calls = 0_8
  INTEGER(kind=8), SAVE :: tangent_idx_n = 0_8
  REAL(kind=8), ALLOCATABLE, SAVE :: tangent_data_ctx(:, :)
  REAL(kind=8), ALLOCATABLE, SAVE :: tangent_cell_ctx(:, :)
  PetscInt, ALLOCATABLE, SAVE :: tangent_idx(:)
  REAL(kind=8), ALLOCATABLE, SAVE :: tangent_work_x(:), tangent_work_y(:)
  PetscScalar, ALLOCATABLE, SAVE :: tangent_work_xvals(:), tangent_work_yvals(:)
  REAL(kind=8), ALLOCATABLE, SAVE :: tangent_work_dw(:)
  LOGICAL, SAVE :: tangent_baseline_ready = .FALSE.
  REAL(kind=8), ALLOCATABLE, SAVE :: baseline_fluxres(:,:), baseline_faceflux(:,:), baseline_cellflux(:,:)
  REAL(kind=8), ALLOCATABLE, SAVE :: baseline_cellstate(:,:), baseline_cellprim(:,:)
  REAL(kind=8), ALLOCATABLE, SAVE :: stack_fluxres(:,:), stack_faceflux(:,:), stack_cellflux(:,:)
  REAL(kind=8), ALLOCATABLE, SAVE :: stack_cellstate(:,:), stack_cellprim(:,:)

CONTAINS

  REAL(kind=8) FUNCTION SAFE_NORM2_2D(arr)
    IMPLICIT NONE
    REAL(kind=8), INTENT(IN) :: arr(:, :)
    SAFE_NORM2_2D = SQRT(SUM(arr*arr))
  END FUNCTION SAFE_NORM2_2D

  SUBROUTINE ENSURE_TANGENT_WORKSPACE(nloc)
    IMPLICIT NONE
    INTEGER(kind=8), INTENT(IN) :: nloc
    IF (.NOT. ALLOCATED(tangent_work_x) .OR. SIZE(tangent_work_x) /= nloc) THEN
      IF (ALLOCATED(tangent_work_x)) DEALLOCATE(tangent_work_x)
      IF (ALLOCATED(tangent_work_y)) DEALLOCATE(tangent_work_y)
      IF (ALLOCATED(tangent_work_xvals)) DEALLOCATE(tangent_work_xvals)
      IF (ALLOCATED(tangent_work_yvals)) DEALLOCATE(tangent_work_yvals)
      IF (ALLOCATED(tangent_work_dw)) DEALLOCATE(tangent_work_dw)
      ALLOCATE(tangent_work_x(nloc), tangent_work_y(nloc), &
     &         tangent_work_xvals(nloc), tangent_work_yvals(nloc), &
     &         tangent_work_dw(nloc))
    END IF
  END SUBROUTINE ENSURE_TANGENT_WORKSPACE

  SUBROUTINE ENSURE_2D_SHAPE(buf, n1, n2)
    IMPLICIT NONE
    REAL(kind=8), ALLOCATABLE, INTENT(INOUT) :: buf(:, :)
    INTEGER, INTENT(IN) :: n1, n2
    IF (.NOT. ALLOCATED(buf) .OR. SIZE(buf,1) /= n1 .OR. SIZE(buf,2) /= n2) THEN
      IF (ALLOCATED(buf)) DEALLOCATE(buf)
      ALLOCATE(buf(n1, n2))
    END IF
  END SUBROUTINE ENSURE_2D_SHAPE

  SUBROUTINE TANGENT_SANDBOX_CAPTURE_BASELINE()
    IMPLICIT NONE
    IF (ALLOCATED(fluxresiduals_slnd)) THEN
      CALL ENSURE_2D_SHAPE(baseline_fluxres, SIZE(fluxresiduals_slnd,1), SIZE(fluxresiduals_slnd,2))
      baseline_fluxres = fluxresiduals_slnd
    END IF
    IF (ALLOCATED(facefluxes_slnd)) THEN
      CALL ENSURE_2D_SHAPE(baseline_faceflux, SIZE(facefluxes_slnd,1), SIZE(facefluxes_slnd,2))
      baseline_faceflux = facefluxes_slnd
    END IF
    IF (ALLOCATED(cellfluxes_slnd)) THEN
      CALL ENSURE_2D_SHAPE(baseline_cellflux, SIZE(cellfluxes_slnd,1), SIZE(cellfluxes_slnd,2))
      baseline_cellflux = cellfluxes_slnd
    END IF
    IF (ALLOCATED(cellstate_slnd)) THEN
      CALL ENSURE_2D_SHAPE(baseline_cellstate, SIZE(cellstate_slnd,1), SIZE(cellstate_slnd,2))
      baseline_cellstate = cellstate_slnd
    END IF
    IF (ALLOCATED(cellprimitives_slnd)) THEN
      CALL ENSURE_2D_SHAPE(baseline_cellprim, SIZE(cellprimitives_slnd,1), SIZE(cellprimitives_slnd,2))
      baseline_cellprim = cellprimitives_slnd
    END IF
    tangent_baseline_ready = .TRUE.
  END SUBROUTINE TANGENT_SANDBOX_CAPTURE_BASELINE

  SUBROUTINE TANGENT_SANDBOX_PUSH()
    IMPLICIT NONE
    IF (ALLOCATED(fluxresiduals_slnd)) THEN
      CALL ENSURE_2D_SHAPE(stack_fluxres, SIZE(fluxresiduals_slnd,1), SIZE(fluxresiduals_slnd,2))
      stack_fluxres = fluxresiduals_slnd
    END IF
    IF (ALLOCATED(facefluxes_slnd)) THEN
      CALL ENSURE_2D_SHAPE(stack_faceflux, SIZE(facefluxes_slnd,1), SIZE(facefluxes_slnd,2))
      stack_faceflux = facefluxes_slnd
    END IF
    IF (ALLOCATED(cellfluxes_slnd)) THEN
      CALL ENSURE_2D_SHAPE(stack_cellflux, SIZE(cellfluxes_slnd,1), SIZE(cellfluxes_slnd,2))
      stack_cellflux = cellfluxes_slnd
    END IF
    IF (ALLOCATED(cellstate_slnd)) THEN
      CALL ENSURE_2D_SHAPE(stack_cellstate, SIZE(cellstate_slnd,1), SIZE(cellstate_slnd,2))
      stack_cellstate = cellstate_slnd
    END IF
    IF (ALLOCATED(cellprimitives_slnd)) THEN
      CALL ENSURE_2D_SHAPE(stack_cellprim, SIZE(cellprimitives_slnd,1), SIZE(cellprimitives_slnd,2))
      stack_cellprim = cellprimitives_slnd
    END IF
  END SUBROUTINE TANGENT_SANDBOX_PUSH

  SUBROUTINE TANGENT_SANDBOX_APPLY_BASELINE()
    IMPLICIT NONE
    IF (.NOT. tangent_baseline_ready) RETURN
    IF (ALLOCATED(fluxresiduals_slnd) .AND. ALLOCATED(baseline_fluxres)) THEN
      IF (SIZE(fluxresiduals_slnd,1) == SIZE(baseline_fluxres,1) .AND. SIZE(fluxresiduals_slnd,2) == SIZE(baseline_fluxres,2)) THEN
        fluxresiduals_slnd = baseline_fluxres
      END IF
    END IF
    IF (ALLOCATED(facefluxes_slnd) .AND. ALLOCATED(baseline_faceflux)) THEN
      IF (SIZE(facefluxes_slnd,1) == SIZE(baseline_faceflux,1) .AND. SIZE(facefluxes_slnd,2) == SIZE(baseline_faceflux,2)) THEN
        facefluxes_slnd = baseline_faceflux
      END IF
    END IF
    IF (ALLOCATED(cellfluxes_slnd) .AND. ALLOCATED(baseline_cellflux)) THEN
      IF (SIZE(cellfluxes_slnd,1) == SIZE(baseline_cellflux,1) .AND. SIZE(cellfluxes_slnd,2) == SIZE(baseline_cellflux,2)) THEN
        cellfluxes_slnd = baseline_cellflux
      END IF
    END IF
    IF (ALLOCATED(cellstate_slnd) .AND. ALLOCATED(baseline_cellstate)) THEN
      IF (SIZE(cellstate_slnd,1) == SIZE(baseline_cellstate,1) .AND. SIZE(cellstate_slnd,2) == SIZE(baseline_cellstate,2)) THEN
        cellstate_slnd = baseline_cellstate
      END IF
    END IF
    IF (ALLOCATED(cellprimitives_slnd) .AND. ALLOCATED(baseline_cellprim)) THEN
      IF (SIZE(cellprimitives_slnd,1) == SIZE(baseline_cellprim,1) .AND. &
     &    SIZE(cellprimitives_slnd,2) == SIZE(baseline_cellprim,2)) THEN
        cellprimitives_slnd = baseline_cellprim
      END IF
    END IF
  END SUBROUTINE TANGENT_SANDBOX_APPLY_BASELINE

  SUBROUTINE TANGENT_SANDBOX_POP()
    IMPLICIT NONE
    IF (ALLOCATED(fluxresiduals_slnd) .AND. ALLOCATED(stack_fluxres)) THEN
      IF (SIZE(fluxresiduals_slnd,1) == SIZE(stack_fluxres,1) .AND. SIZE(fluxresiduals_slnd,2) == SIZE(stack_fluxres,2)) THEN
        fluxresiduals_slnd = stack_fluxres
      END IF
    END IF
    IF (ALLOCATED(facefluxes_slnd) .AND. ALLOCATED(stack_faceflux)) THEN
      IF (SIZE(facefluxes_slnd,1) == SIZE(stack_faceflux,1) .AND. SIZE(facefluxes_slnd,2) == SIZE(stack_faceflux,2)) THEN
        facefluxes_slnd = stack_faceflux
      END IF
    END IF
    IF (ALLOCATED(cellfluxes_slnd) .AND. ALLOCATED(stack_cellflux)) THEN
      IF (SIZE(cellfluxes_slnd,1) == SIZE(stack_cellflux,1) .AND. SIZE(cellfluxes_slnd,2) == SIZE(stack_cellflux,2)) THEN
        cellfluxes_slnd = stack_cellflux
      END IF
    END IF
    IF (ALLOCATED(cellstate_slnd) .AND. ALLOCATED(stack_cellstate)) THEN
      IF (SIZE(cellstate_slnd,1) == SIZE(stack_cellstate,1) .AND. SIZE(cellstate_slnd,2) == SIZE(stack_cellstate,2)) THEN
        cellstate_slnd = stack_cellstate
      END IF
    END IF
    IF (ALLOCATED(cellprimitives_slnd) .AND. ALLOCATED(stack_cellprim)) THEN
      IF (SIZE(cellprimitives_slnd,1) == SIZE(stack_cellprim,1) .AND. &
     &    SIZE(cellprimitives_slnd,2) == SIZE(stack_cellprim,2)) THEN
        cellprimitives_slnd = stack_cellprim
      END IF
    END IF
  END SUBROUTINE TANGENT_SANDBOX_POP

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

  SUBROUTINE BUILD_PMAT_FROM_PC_BLOCKS(nloc, Pmat, pmat_ready, ierr)
    IMPLICIT NONE
    INTEGER(kind=8), INTENT(IN) :: nloc
    Mat, INTENT(OUT) :: Pmat
    LOGICAL, INTENT(OUT) :: pmat_ready
    PetscErrorCode, INTENT(OUT) :: ierr
    INTEGER(kind=8) :: ncells, c, p, j, ir, jc, nblkrow
    PetscInt :: n_petsc
    PetscInt, ALLOCATABLE :: d_nnz(:), row_idx(:), col_idx(:)
    PetscScalar :: row_vals(5)
    LOGICAL :: use_a0

    ierr = 0
    pmat_ready = .FALSE.

    IF (.NOT. ALLOCATED(pc_row_ptr) .OR. .NOT. ALLOCATED(pc_col_ind) .OR. .NOT. ALLOCATED(pc_diag_pos)) RETURN
    use_a0 = pc_a0_ready .AND. ALLOCATED(pc_a0_blk)
    IF ((.NOT. use_a0) .AND. (.NOT. ALLOCATED(pc_blk))) RETURN

    ncells = SIZE(pc_diag_pos, kind=8)
    IF (ncells <= 0_8) RETURN
    IF (nloc /= 5_8 * ncells) RETURN
    IF (SIZE(pc_row_ptr, kind=8) < ncells + 1_8) RETURN

    n_petsc = INT(nloc, KIND=KIND(n_petsc))
    ALLOCATE(d_nnz(n_petsc))
    d_nnz = 1
    DO c = 1_8, ncells
      nblkrow = MAX(1_8, pc_row_ptr(c+1_8) - pc_row_ptr(c))
      DO ir = 1_8, 5_8
        d_nnz(INT((c-1_8)*5_8 + ir, KIND=KIND(d_nnz(1)))) = INT(5_8*nblkrow, KIND=KIND(d_nnz(1)))
      END DO
    END DO

    CALL MatCreateSeqAIJ(PETSC_COMM_SELF, n_petsc, n_petsc, 0, d_nnz, Pmat, ierr)
    DEALLOCATE(d_nnz)
    IF (ierr /= 0) RETURN

    ALLOCATE(row_idx(1), col_idx(5))
    DO c = 1_8, ncells
      DO p = pc_row_ptr(c), pc_row_ptr(c+1_8)-1_8
        j = pc_col_ind(p)
        IF (j < 1_8 .OR. j > ncells) CYCLE

        DO jc = 1_8, 5_8
          col_idx(jc) = INT((j-1_8)*5_8 + (jc-1_8), KIND=KIND(col_idx(jc)))
        END DO
        DO ir = 1_8, 5_8
          row_idx(1) = INT((c-1_8)*5_8 + (ir-1_8), KIND=KIND(row_idx(1)))
          DO jc = 1_8, 5_8
            IF (use_a0) THEN
              row_vals(jc) = pc_a0_blk(p, ir, jc)
            ELSE
              row_vals(jc) = pc_blk(p, ir, jc)
            END IF
          END DO
          CALL MatSetValues(Pmat, 1, row_idx, 5, col_idx, row_vals, INSERT_VALUES, ierr)
          IF (ierr /= 0) THEN
            DEALLOCATE(row_idx, col_idx)
            RETURN
          END IF
        END DO
      END DO
    END DO
    DEALLOCATE(row_idx, col_idx)

    CALL MatAssemblyBegin(Pmat, MAT_FINAL_ASSEMBLY, ierr)
    IF (ierr /= 0) RETURN
    CALL MatAssemblyEnd(Pmat, MAT_FINAL_ASSEMBLY, ierr)
    IF (ierr /= 0) RETURN

    pmat_ready = .TRUE.
  END SUBROUTINE BUILD_PMAT_FROM_PC_BLOCKS

  SUBROUTINE TANGENT_PETSC_MATMULT_IMPL(A, X, Y, ierr)
    IMPLICIT NONE
    Mat :: A, Pmat
    Vec :: X, Y
    PetscErrorCode :: ierr
    INTEGER(kind=8) :: nloc, i
    LOGICAL :: had_bad
    REAL(kind=8) :: xin_norm, yout_norm, flux_pre, flux_post, point_pre, point_post
    REAL(kind=8) :: max_abs_raw

    ierr = 0
    nloc = tangent_n_ctx
    IF (nloc <= 0_8) RETURN

    IF (tangent_idx_n /= nloc .OR. .NOT. ALLOCATED(tangent_idx)) THEN
      ierr = 1
      RETURN
    END IF
    CALL ENSURE_TANGENT_WORKSPACE(nloc)
    tangent_work_x = 0.0_8
    tangent_work_y = 0.0_8
    tangent_work_xvals = 0.0_8
    tangent_work_yvals = 0.0_8

    CALL VecGetValues(X, INT(nloc, kind=KIND(1)), tangent_idx, tangent_work_xvals, ierr)
    IF (ierr /= 0) THEN
      RETURN
    END IF
    tangent_work_x = tangent_work_xvals

    had_bad = .FALSE.
    DO i = 1_8, nloc
      IF (tangent_work_x(i) /= tangent_work_x(i) .OR. ABS(tangent_work_x(i)) > 1.0d120) THEN
        tangent_work_x(i) = 0.0_8
        had_bad = .TRUE.
      END IF
    END DO

    flux_pre = -1.0_8
    flux_post = -1.0_8
    point_pre = -1.0_8
    point_post = -1.0_8
    IF (ALLOCATED(fluxresiduals_slnd)) flux_pre = SAFE_NORM2_2D(fluxresiduals_slnd)
    IF (ALLOCATED(point_updated)) point_pre = SAFE_NORM2_2D(point_updated)
    CALL TANGENT_SANDBOX_PUSH()
    CALL TANGENT_SANDBOX_APPLY_BASELINE()
    CALL APPLY_TANGENT_OPERATOR_WORK(tangent_data_ctx, tangent_cell_ctx, nloc, &
   &     tangent_work_x, 0.0_8, tangent_work_y, tangent_work_dw)
    CALL APPLY_RESIDUAL_SCALING(tangent_work_y)
    CALL TANGENT_SANDBOX_POP()
    IF (ALLOCATED(fluxresiduals_slnd)) flux_post = SAFE_NORM2_2D(fluxresiduals_slnd)
    IF (ALLOCATED(point_updated)) point_post = SAFE_NORM2_2D(point_updated)
    IF (tangent_matmult_calls == 0_8) THEN
      PRINT *, '[PETSC-MATMULT-DBG] raw yout(1:5)=', tangent_work_y(1:MIN(5_8, nloc))
    END IF

    max_abs_raw = 0.0_8
    DO i = 1_8, nloc
      max_abs_raw = MAX(max_abs_raw, ABS(tangent_work_y(i)))
      IF (tangent_work_y(i) /= tangent_work_y(i) .OR. ABS(tangent_work_y(i)) > 1.0d120) THEN
        IF (.NOT. had_bad .AND. tangent_matmult_calls <= 2_8) THEN
          PRINT *, '[PETSC-MATMULT-DBG] first bad yout idx=', i, ' value=', tangent_work_y(i)
        END IF
        tangent_work_y(i) = 0.0_8
        had_bad = .TRUE.
      END IF
    END DO

    tangent_work_yvals = tangent_work_y
    CALL VecSetValues(Y, INT(nloc, kind=KIND(1)), tangent_idx, tangent_work_yvals, INSERT_VALUES, ierr)
    CALL VecAssemblyBegin(Y, ierr)
    CALL VecAssemblyEnd(Y, ierr)

    tangent_matmult_calls = tangent_matmult_calls + 1_8
    IF (tangent_matmult_calls <= 3_8) THEN
      xin_norm = SQRT(SUM(tangent_work_x*tangent_work_x))
      yout_norm = SQRT(SUM(tangent_work_y*tangent_work_y))
      PRINT *, '[PETSC-MATMULT] call=', tangent_matmult_calls, ' ||x||=', xin_norm, ' ||Ax||=', yout_norm, ' sanitized=', had_bad
      PRINT *, '[PETSC-SNAP] flux pre/post=', flux_pre, flux_post, &
     &         ' point pre/post=', point_pre, point_post, &
     &         ' max|Ax_raw|=', max_abs_raw
    END IF
  END SUBROUTINE TANGENT_PETSC_MATMULT_IMPL

  SUBROUTINE SOLVE_TANGENT_PETSC_GMRES(data_4d137, cellprimitives, n, b, x, tol, maxiter, m_restart, info)
    IMPLICIT NONE
    REAL(kind=8), INTENT(IN) :: data_4d137(1, 4), cellprimitives(:, :), b(n), tol
    REAL(kind=8), INTENT(INOUT) :: x(n)
    INTEGER(kind=8), INTENT(IN) :: n, maxiter, m_restart
    INTEGER, INTENT(OUT) :: info

    Mat :: A, Pmat
    Vec :: vb, vx, vtmp, vr
    KSP :: ksp
    PC :: pc
    KSPConvergedReason :: reason
    PetscErrorCode :: ierr
    PetscInt :: n_petsc, max_it_petsc, restart_petsc, its, i
    PetscReal :: rnorm, rnorm_true, bnorm, xnorm0, axnorm0, axnormb, az_pre1, az_pre2, az_pre3, az_post
    PetscInt, ALLOCATABLE :: idx(:)
    PetscScalar, ALLOCATABLE :: vals(:)
    REAL(kind=8), ALLOCATABLE :: zprobe(:), az(:), b_work(:)
    EXTERNAL TANGENT_PETSC_MATMULT
    LOGICAL :: pmat_ready, pmat_created

    info = 2
    pmat_created = .FALSE.
    IF (n <= 0_8) RETURN

    CALL TANGENT_PETSC_SET_CONTEXT(data_4d137, cellprimitives, n)
    ALLOCATE(b_work(n))
    b_work = b
    CALL APPLY_RESIDUAL_SCALING(b_work)
    ALLOCATE(zprobe(n), az(n))
    zprobe = 0.0_8
    CALL APPLY_TANGENT_OPERATOR_WORK(tangent_data_ctx, tangent_cell_ctx, n, zprobe, 0.0_8, az, tangent_work_dw)
    CALL APPLY_RESIDUAL_SCALING(az)
    az_pre1 = SQRT(SUM(az*az))
    PRINT *, '[PETSC-PREINIT-1] ||A*0||=', az_pre1
    CALL APPLY_TANGENT_OPERATOR_WORK(tangent_data_ctx, tangent_cell_ctx, n, zprobe, 0.0_8, az, tangent_work_dw)
    CALL APPLY_RESIDUAL_SCALING(az)
    az_pre2 = SQRT(SUM(az*az))
    PRINT *, '[PETSC-PREINIT-2] ||A*0||=', az_pre2
    CALL APPLY_TANGENT_OPERATOR_WORK(tangent_data_ctx, tangent_cell_ctx, n, zprobe, 0.0_8, az, tangent_work_dw)
    CALL APPLY_RESIDUAL_SCALING(az)
    az_pre3 = SQRT(SUM(az*az))
    PRINT *, '[PETSC-PREINIT-3] ||A*0||=', az_pre3
    CALL TANGENT_SANDBOX_CAPTURE_BASELINE()
    PRINT *, '[PETSC-SANDBOX] baseline captured'

    IF (.NOT. petsc_initialized_local) THEN
      CALL PetscInitialize(ierr)
      IF (ierr /= 0) THEN
        PRINT *, '[PETSC] initialize failed, ierr=', ierr
        RETURN
      END IF
      petsc_initialized_local = .TRUE.
    END IF
    CALL APPLY_TANGENT_OPERATOR_WORK(tangent_data_ctx, tangent_cell_ctx, n, zprobe, 0.0_8, az, tangent_work_dw)
    CALL APPLY_RESIDUAL_SCALING(az)
    az_post = SQRT(SUM(az*az))
    PRINT *, '[PETSC-POSTINIT] ||A*0||=', az_post

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
    CALL BUILD_PMAT_FROM_PC_BLOCKS(n, Pmat, pmat_ready, ierr)
    IF (ierr == 0 .AND. pmat_ready) THEN
      pmat_created = .TRUE.
      PRINT *, '[PETSC-PMAT] using explicit AIJ Pmat from PC blocks'
    ELSE
      pmat_ready = .FALSE.
      ierr = 0
      PRINT *, '[PETSC-PMAT] fallback: no explicit Pmat available'
    END IF

    CALL VecCreateSeq(PETSC_COMM_SELF, n_petsc, vb, ierr)
    CALL VecDuplicate(vb, vx, ierr)

    ALLOCATE(idx(n_petsc), vals(n_petsc))
    DO i = 1, n_petsc
      idx(i) = i - 1
    END DO
    vals = b_work(1:n)
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
    IF (pmat_ready) THEN
      CALL KSPSetOperators(ksp, A, Pmat, ierr)
    ELSE
      CALL KSPSetOperators(ksp, A, A, ierr)
    END IF
    ! Keep a safe default; actual KSP/PC is overridden from runtime options.
    CALL KSPSetType(ksp, KSPGMRES, ierr)
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

    CALL VecDuplicate(vb, vr, ierr)
    CALL MatMult(A, vx, vr, ierr)
    CALL VecAYPX(vr, -1.0d0, vb, ierr) ! vr = b - A*x
    CALL VecNorm(vr, NORM_2, rnorm_true, ierr)
    CALL VecDestroy(vr, ierr)

    bnorm = SQRT(SUM(b_work*b_work))
    PRINT *, '[PETSC-GMRES] reason=', reason, ' its=', its, ' ksp||r||/||b||=', &
   & rnorm / MAX(1.0d-30, bnorm), ' true||r||/||b||=', rnorm_true / MAX(1.0d-30, bnorm)

    IF (reason > 0) THEN
      info = 0
    ELSE
      info = 1
    END IF

    CALL KSPDestroy(ksp, ierr)
    CALL VecDestroy(vb, ierr)
    CALL VecDestroy(vx, ierr)
    IF (pmat_created) CALL MatDestroy(Pmat, ierr)
    CALL MatDestroy(A, ierr)
    DEALLOCATE(b_work)
    DEALLOCATE(zprobe, az)
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
    IF (ALLOCATED(tangent_work_x)) DEALLOCATE(tangent_work_x)
    IF (ALLOCATED(tangent_work_y)) DEALLOCATE(tangent_work_y)
    IF (ALLOCATED(tangent_work_xvals)) DEALLOCATE(tangent_work_xvals)
    IF (ALLOCATED(tangent_work_yvals)) DEALLOCATE(tangent_work_yvals)
    IF (ALLOCATED(tangent_work_dw)) DEALLOCATE(tangent_work_dw)
    IF (ALLOCATED(baseline_fluxres)) DEALLOCATE(baseline_fluxres)
    IF (ALLOCATED(baseline_faceflux)) DEALLOCATE(baseline_faceflux)
    IF (ALLOCATED(baseline_cellflux)) DEALLOCATE(baseline_cellflux)
    IF (ALLOCATED(baseline_cellstate)) DEALLOCATE(baseline_cellstate)
    IF (ALLOCATED(baseline_cellprim)) DEALLOCATE(baseline_cellprim)
    IF (ALLOCATED(stack_fluxres)) DEALLOCATE(stack_fluxres)
    IF (ALLOCATED(stack_faceflux)) DEALLOCATE(stack_faceflux)
    IF (ALLOCATED(stack_cellflux)) DEALLOCATE(stack_cellflux)
    IF (ALLOCATED(stack_cellstate)) DEALLOCATE(stack_cellstate)
    IF (ALLOCATED(stack_cellprim)) DEALLOCATE(stack_cellprim)
    tangent_idx_n = 0_8
    tangent_n_ctx = 0_8
    tangent_baseline_ready = .FALSE.
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
