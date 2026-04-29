!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief  Functionality for atom centered symmetry functions
!>         for neural network potentials
!> \author Christoph Schran (christoph.schran@rub.de)
!> \date   2020-10-10
! **************************************************************************************************
MODULE nnp_acsf
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE mathconstants,                   ONLY: pi
   USE message_passing,                 ONLY: mp_para_env_type
   USE nnp_environment_types,           ONLY: nnp_acsf_ang_type,&
                                              nnp_acsf_rad_type,&
                                              nnp_cut_cos,&
                                              nnp_cut_tanh,&
                                              nnp_env_get,&
                                              nnp_neighbor_type,&
                                              nnp_type
   USE periodic_table,                  ONLY: get_ptable_info
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'nnp_acsf'
   ! Public subroutines ***
   PUBLIC :: nnp_calc_acsf, &
             nnp_init_acsf_groups, &
             nnp_prepare_cell_list_cache, &
             nnp_sort_acsf, &
             nnp_sort_ele, &
             nnp_write_acsf

! **************************************************************************************************
!> \brief Module-level cache for the pre-expanded-image cell list used by the O(N) neighbour query.
!>        Each (atom j, integer PBC shift (tx,ty,tz)) becomes one image entry; image positions are
!>        stored explicitly so the inner query loop does no pbc() calls.  Grid is coarse
!>        (bin_width ~ list_cutoff, bin_span = 1) so each centre atom walks 27 bins.  Rebuild is
!>        triggered only by Verlet-skin violation (any atom moved > 0.5*skin) or by changes to
!>        num_atoms / cell vectors / cutoff.
!> \author Claudio Malvino
!> \date   2026-04-24
! **************************************************************************************************
   TYPE, PRIVATE :: nnp_cell_list_cache_type
      LOGICAL                              :: initialized   = .FALSE.
      INTEGER                              :: num_atoms     = -1
      INTEGER                              :: n_images      = 0
      INTEGER                              :: n_cells       = 1
      INTEGER, DIMENSION(3)                :: image_copies  = 0
      INTEGER, DIMENSION(3)                :: nbin          = 1
      INTEGER, DIMENSION(3)                :: bin_span      = 0
      INTEGER, DIMENSION(3)                :: perd          = 0
      LOGICAL                              :: orthorhombic  = .FALSE.
      REAL(KIND=dp)                        :: exact_cutoff  = -1.0_dp
      REAL(KIND=dp)                        :: list_cutoff   = -1.0_dp
      REAL(KIND=dp)                        :: verlet_skin   = 0.0_dp
      REAL(KIND=dp), DIMENSION(3)          :: lower         = 0.0_dp
      REAL(KIND=dp), DIMENSION(3)          :: upper         = 0.0_dp
      REAL(KIND=dp), DIMENSION(3)          :: bin_width     = 0.0_dp
      REAL(KIND=dp), DIMENSION(3, 3)       :: hmat          = 0.0_dp
      REAL(KIND=dp), DIMENSION(3, 3)       :: h_inv         = 0.0_dp
      ! Primary (PBC-wrapped) atom positions — refreshed every prepare call
      REAL(KIND=dp), ALLOCATABLE           :: coord_primary(:, :)         ! (3, num_atoms)
      ! Snapshot of coord_primary at the last rebuild (Verlet reference)
      REAL(KIND=dp), ALLOCATABLE           :: ref_coord_primary(:, :)     ! (3, num_atoms)
      ! Image table: one entry per (atom j, integer PBC shift) pair
      INTEGER, ALLOCATABLE                 :: image_atom(:)               ! (n_images)
      INTEGER, ALLOCATABLE                 :: image_shift(:, :)           ! (3, n_images)
      REAL(KIND=dp), ALLOCATABLE           :: image_translation(:, :)     ! (3, n_images)
      ! Linked-list cell grid: head(n_cells) → first image in bin, next(n_images) chains images
      INTEGER, ALLOCATABLE                 :: head(:)                     ! (n_cells)
      INTEGER, ALLOCATABLE                 :: next(:)                     ! (n_images)
   END TYPE nnp_cell_list_cache_type

! **************************************************************************************************
!> \brief Module-level cache for Hermite Cubic Splines to replace expensive intrinsic functions.
!>        Stores coefficients (a, b, c, d) for Horner's method evaluation.
! **************************************************************************************************
   TYPE, PRIVATE :: nnp_spline_cache_type
      LOGICAL                                :: initialized = .FALSE.
      INTEGER                                :: n_points = 5000
      REAL(KIND=dp)                          :: dx_exp, inv_dx_exp, min_exp
      REAL(KIND=dp)                          :: dx_cos, inv_dx_cos, min_cos, max_cos
      REAL(KIND=dp)                          :: dx_tanh, inv_dx_tanh, min_tanh, max_tanh
      
      ! Dimension 1 is contiguous in memory (a, b, c, d). Dimension 2 is the grid index.
      REAL(KIND=dp), ALLOCATABLE             :: coef_exp(:,:)
      REAL(KIND=dp), ALLOCATABLE             :: coef_cos(:,:)
      REAL(KIND=dp), ALLOCATABLE             :: coef_tanh(:,:)
   END TYPE nnp_spline_cache_type

   TYPE(nnp_spline_cache_type), SAVE, PRIVATE :: spline_cache
   TYPE(nnp_cell_list_cache_type), SAVE, PRIVATE :: cell_list_cache

   ! Persistent per-centre neighbour buffer, reused across all atoms in a force eval.
   ! Without this each centre atom would do 9 ALLOCATEs sized in num_atoms — at ~6k atoms those
   ! cross glibc's mmap threshold and become syscalls, giving O(N^2) malloc traffic per step.
   ! Sized to the MAX (n_symfgrp_rad / n_symfgrp_ang) over element types so the same buffer
   ! works for every centre regardless of which element it is.
   TYPE(nnp_neighbor_type), SAVE, PRIVATE, TARGET :: persistent_neighbor
   INTEGER, SAVE, PRIVATE                          :: persistent_n_capacity = -1
   INTEGER, SAVE, PRIVATE                          :: persistent_n_rad_grp = -1
   INTEGER, SAVE, PRIVATE                          :: persistent_n_ang_grp = -1
   ! Module-level scratch buffers to avoid repeated allocate/deallocate in hot loops
   INTEGER, SAVE, PRIVATE                          :: scratch_max_sym = -1
   REAL(KIND=dp), ALLOCATABLE, SAVE, TARGET, PRIVATE :: scratch_sym(:)
   REAL(KIND=dp), ALLOCATABLE, SAVE, TARGET, PRIVATE :: scratch_forcetmp(:,:)
   REAL(KIND=dp), ALLOCATABLE, SAVE, TARGET, PRIVATE :: scratch_force3tmp(:,:,:)

CONTAINS

! **************************************************************************************************
!> \brief Calculate atom centered symmetry functions for given atom i
!> \param nnp ...
!> \param i ...
!> \param dsymdxyz ...
!> \param stress ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_calc_acsf(nnp, i, dsymdxyz, stress, active_atoms, active_flag, n_active)

      TYPE(nnp_type), INTENT(INOUT), POINTER             :: nnp
      INTEGER, INTENT(IN)                                :: i
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT), &
         OPTIONAL                                        :: dsymdxyz, stress
      INTEGER, DIMENSION(:), OPTIONAL, INTENT(INOUT)     :: active_atoms, active_flag
      INTEGER, OPTIONAL, INTENT(OUT)                     :: n_active

      CHARACTER(len=*), PARAMETER                        :: routineN = 'nnp_calc_acsf'

      INTEGER                                            :: handle, handle_sf, ind, j, jj, k, kk, l, &
                                                            m, off, s, sf
      REAL(KIND=dp)                                      :: r1, r2, r3
      REAL(KIND=dp), DIMENSION(:), POINTER               :: symtmp
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: forcetmp
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: force3tmp
      REAL(KIND=dp), DIMENSION(3)                        :: rvect1, rvect2, rvect3
      TYPE(nnp_neighbor_type), POINTER                   :: neighbor

      CALL timeset(routineN, handle)

      !determine index of atom type
      ind = nnp%ele_ind(i)

      ! Reuse the module-level persistent neighbour buffer (sized once per force eval by
      ! nnp_neighbor_ensure_persistent).  Only the counts are reset between centre atoms;
      ! no allocations happen in the per-atom hot path.
      neighbor => persistent_neighbor
      neighbor%n_rad(:) = 0
      neighbor%n_ang1(:) = 0
      neighbor%n_ang2(:) = 0
      CALL nnp_compute_neighbors(nnp, neighbor, i)

      ! Build active-atoms list: centre i + every distinct neighbour touched by the rad/ang
      ! lists.  active_flag is the dedup map (caller resets touched entries to 0 between centres).
      ! Restricts the O(num_atoms) dsymdxyz / force-accumulation loops to the actual neighbour
      ! set — drops O(N^2) per force eval to O(N) at fixed density.
      IF (PRESENT(active_atoms) .AND. PRESENT(active_flag) .AND. PRESENT(n_active)) THEN
         n_active = 1
         active_atoms(1) = i
         active_flag(i) = 1
         DO s = 1, nnp%rad(ind)%n_symfgrp
            DO j = 1, neighbor%n_rad(s)
               jj = neighbor%ind_rad(j, s)
               IF (active_flag(jj) == 0) THEN
                  n_active = n_active + 1
                  active_atoms(n_active) = jj
                  active_flag(jj) = 1
               END IF
            END DO
         END DO
         DO s = 1, nnp%ang(ind)%n_symfgrp
            DO j = 1, neighbor%n_ang1(s)
               jj = neighbor%ind_ang1(j, s)
               IF (active_flag(jj) == 0) THEN
                  n_active = n_active + 1
                  active_atoms(n_active) = jj
                  active_flag(jj) = 1
               END IF
            END DO
            DO j = 1, neighbor%n_ang2(s)
               jj = neighbor%ind_ang2(j, s)
               IF (active_flag(jj) == 0) THEN
                  n_active = n_active + 1
                  active_atoms(n_active) = jj
                  active_flag(jj) = 1
               END IF
            END DO
         END DO
      END IF

      ! Reset y:
      nnp%rad(ind)%y = 0.0_dp
      nnp%ang(ind)%y = 0.0_dp

      !calc forces
      IF (PRESENT(dsymdxyz)) THEN
         !loop over radial sym fnct grps
         CALL timeset('nnp_acsf_radial', handle_sf)
         DO s = 1, nnp%rad(ind)%n_symfgrp
            ! reuse module-level scratch buffers to avoid repeated alloc/dealloc
            IF (.NOT. ALLOCATED(scratch_sym)) THEN
               ALLOCATE(scratch_sym(nnp%rad(ind)%symfgrp(s)%n_symf))
               ALLOCATE(scratch_forcetmp(3, nnp%rad(ind)%symfgrp(s)%n_symf))
               ALLOCATE(scratch_force3tmp(3,3,nnp%rad(ind)%symfgrp(s)%n_symf))
               scratch_max_sym = SIZE(scratch_sym)
            ELSE IF (scratch_max_sym < nnp%rad(ind)%symfgrp(s)%n_symf) THEN
               ! grow to fit if necessary
               DEALLOCATE(scratch_sym); DEALLOCATE(scratch_forcetmp); DEALLOCATE(scratch_force3tmp)
               ALLOCATE(scratch_sym(nnp%rad(ind)%symfgrp(s)%n_symf))
               ALLOCATE(scratch_forcetmp(3, nnp%rad(ind)%symfgrp(s)%n_symf))
               ALLOCATE(scratch_force3tmp(3,3,nnp%rad(ind)%symfgrp(s)%n_symf))
               scratch_max_sym = SIZE(scratch_sym)
            END IF
            symtmp => scratch_sym(1:nnp%rad(ind)%symfgrp(s)%n_symf)
            forcetmp => scratch_forcetmp(:, 1:nnp%rad(ind)%symfgrp(s)%n_symf)
            !loop over associated neighbors
            DO j = 1, neighbor%n_rad(s)
               rvect1 = neighbor%dist_rad(1:3, j, s)
               r1 = neighbor%dist_rad(4, j, s)
               CALL nnp_calc_rad(nnp, ind, s, rvect1, r1, symtmp, forcetmp)
               jj = neighbor%ind_rad(j, s)
               DO sf = 1, nnp%rad(ind)%symfgrp(s)%n_symf
                  m = nnp%rad(ind)%symfgrp(s)%symf(sf)
                  ! update forces into dsymdxyz
                  DO l = 1, 3
                     dsymdxyz(l, m, i) = dsymdxyz(l, m, i) + forcetmp(l, sf)
                     dsymdxyz(l, m, jj) = dsymdxyz(l, m, jj) - forcetmp(l, sf)
                  END DO
                  IF (PRESENT(stress)) THEN
                     DO l = 1, 3
                        stress(:, l, m) = stress(:, l, m) + rvect1(:)*forcetmp(l, sf)
                     END DO
                  END IF
                  nnp%rad(ind)%y(m) = nnp%rad(ind)%y(m) + symtmp(sf)
               END DO
            END DO
            NULLIFY(symtmp)
            NULLIFY(forcetmp)
         END DO
         CALL timestop(handle_sf)
         !loop over angular sym fnct grps
         CALL timeset('nnp_acsf_angular', handle_sf)
         off = nnp%n_rad(ind)
         DO s = 1, nnp%ang(ind)%n_symfgrp
            ! reuse module-level scratch buffers to avoid repeated alloc/dealloc
            IF (.NOT. ALLOCATED(scratch_sym)) THEN
               ALLOCATE(scratch_sym(nnp%ang(ind)%symfgrp(s)%n_symf))
               ALLOCATE(scratch_forcetmp(3, nnp%ang(ind)%symfgrp(s)%n_symf))
               ALLOCATE(scratch_force3tmp(3,3,nnp%ang(ind)%symfgrp(s)%n_symf))
               scratch_max_sym = SIZE(scratch_sym)
            ELSE IF (scratch_max_sym < nnp%ang(ind)%symfgrp(s)%n_symf) THEN
               DEALLOCATE(scratch_sym); DEALLOCATE(scratch_forcetmp); DEALLOCATE(scratch_force3tmp)
               ALLOCATE(scratch_sym(nnp%ang(ind)%symfgrp(s)%n_symf))
               ALLOCATE(scratch_forcetmp(3, nnp%ang(ind)%symfgrp(s)%n_symf))
               ALLOCATE(scratch_force3tmp(3,3,nnp%ang(ind)%symfgrp(s)%n_symf))
               scratch_max_sym = SIZE(scratch_sym)
            END IF
            symtmp => scratch_sym(1:nnp%ang(ind)%symfgrp(s)%n_symf)
            force3tmp => scratch_force3tmp(:,:, 1:nnp%ang(ind)%symfgrp(s)%n_symf)
            !loop over associated neighbors
            IF (nnp%ang(ind)%symfgrp(s)%ele(1) == nnp%ang(ind)%symfgrp(s)%ele(2)) THEN
               DO j = 1, neighbor%n_ang1(s)
                  rvect1 = neighbor%dist_ang1(1:3, j, s)
                  r1 = neighbor%dist_ang1(4, j, s)
                  DO k = j + 1, neighbor%n_ang1(s)
                     rvect2 = neighbor%dist_ang1(1:3, k, s)
                     r2 = neighbor%dist_ang1(4, k, s)
                     rvect3 = rvect2 - rvect1
                     r3 = NORM2(rvect3(:))
                     IF (r3 < nnp%ang(ind)%symfgrp(s)%cutoff) THEN
                        jj = neighbor%ind_ang1(j, s)
                        kk = neighbor%ind_ang1(k, s)
                        CALL nnp_calc_ang(nnp, ind, s, rvect1, rvect2, rvect3, &
                                          r1, r2, r3, symtmp, force3tmp)
                        DO sf = 1, nnp%ang(ind)%symfgrp(s)%n_symf
                           m = off + nnp%ang(ind)%symfgrp(s)%symf(sf)
                           ! update forces into dsymdxy
                           DO l = 1, 3
                              dsymdxyz(l, m, i) = dsymdxyz(l, m, i) &
                                                  + force3tmp(l, 1, sf)
                              dsymdxyz(l, m, jj) = dsymdxyz(l, m, jj) &
                                                   + force3tmp(l, 2, sf)
                              dsymdxyz(l, m, kk) = dsymdxyz(l, m, kk) &
                                                   + force3tmp(l, 3, sf)
                           END DO
                           IF (PRESENT(stress)) THEN
                              DO l = 1, 3
                                 stress(:, l, m) = stress(:, l, m) - rvect1(:)*force3tmp(l, 2, sf)
                                 stress(:, l, m) = stress(:, l, m) - rvect2(:)*force3tmp(l, 3, sf)
                              END DO
                           END IF
                           nnp%ang(ind)%y(m - off) = nnp%ang(ind)%y(m - off) + symtmp(sf)
                        END DO
                     END IF
                  END DO
               END DO
            ELSE
               DO j = 1, neighbor%n_ang1(s)
                  rvect1 = neighbor%dist_ang1(1:3, j, s)
                  r1 = neighbor%dist_ang1(4, j, s)
                  DO k = 1, neighbor%n_ang2(s)
                     rvect2 = neighbor%dist_ang2(1:3, k, s)
                     r2 = neighbor%dist_ang2(4, k, s)
                     rvect3 = rvect2 - rvect1
                     r3 = NORM2(rvect3(:))
                     IF (r3 < nnp%ang(ind)%symfgrp(s)%cutoff) THEN
                        jj = neighbor%ind_ang1(j, s)
                        kk = neighbor%ind_ang1(k, s)
                        CALL nnp_calc_ang(nnp, ind, s, rvect1, rvect2, rvect3, &
                                          r1, r2, r3, symtmp, force3tmp)
                        !loop over associated sym fncts
                        DO sf = 1, nnp%ang(ind)%symfgrp(s)%n_symf
                           m = off + nnp%ang(ind)%symfgrp(s)%symf(sf)
                           jj = neighbor%ind_ang1(j, s)
                           kk = neighbor%ind_ang2(k, s)
                           ! update forces into dsymdxy
                           DO l = 1, 3
                              dsymdxyz(l, m, i) = dsymdxyz(l, m, i) &
                                                  + force3tmp(l, 1, sf)
                              dsymdxyz(l, m, jj) = dsymdxyz(l, m, jj) &
                                                   + force3tmp(l, 2, sf)
                              dsymdxyz(l, m, kk) = dsymdxyz(l, m, kk) &
                                                   + force3tmp(l, 3, sf)
                           END DO
                           IF (PRESENT(stress)) THEN
                              DO l = 1, 3
                                 stress(:, l, m) = stress(:, l, m) - rvect1(:)*force3tmp(l, 2, sf)
                                 stress(:, l, m) = stress(:, l, m) - rvect2(:)*force3tmp(l, 3, sf)
                              END DO
                           END IF
                           nnp%ang(ind)%y(m - off) = nnp%ang(ind)%y(m - off) + symtmp(sf)
                        END DO
                     END IF
                  END DO
               END DO
            END IF
            NULLIFY(symtmp)
            NULLIFY(force3tmp)
         END DO
         CALL timestop(handle_sf)
         !no forces
      ELSE
         !loop over radial sym fnct grps
         CALL timeset('nnp_acsf_radial', handle_sf)
         DO s = 1, nnp%rad(ind)%n_symfgrp
            IF (.NOT. ALLOCATED(scratch_sym)) THEN
               ALLOCATE(scratch_sym(nnp%rad(ind)%symfgrp(s)%n_symf))
               scratch_max_sym = SIZE(scratch_sym)
            ELSE IF (scratch_max_sym < nnp%rad(ind)%symfgrp(s)%n_symf) THEN
               DEALLOCATE(scratch_sym)
               ALLOCATE(scratch_sym(nnp%rad(ind)%symfgrp(s)%n_symf))
               scratch_max_sym = SIZE(scratch_sym)
            END IF
            symtmp => scratch_sym(1:nnp%rad(ind)%symfgrp(s)%n_symf)
            !loop over associated neighbors
            DO j = 1, neighbor%n_rad(s)
               rvect1 = neighbor%dist_rad(1:3, j, s)
               r1 = neighbor%dist_rad(4, j, s)
               CALL nnp_calc_rad(nnp, ind, s, rvect1, r1, symtmp)
               DO sf = 1, nnp%rad(ind)%symfgrp(s)%n_symf
                  m = nnp%rad(ind)%symfgrp(s)%symf(sf)
                  nnp%rad(ind)%y(m) = nnp%rad(ind)%y(m) + symtmp(sf)
               END DO
            END DO
            NULLIFY(symtmp)
         END DO
         CALL timestop(handle_sf)
         !loop over angular sym fnct grps
         CALL timeset('nnp_acsf_angular', handle_sf)
         off = nnp%n_rad(ind)
         DO s = 1, nnp%ang(ind)%n_symfgrp
            IF (.NOT. ALLOCATED(scratch_sym)) THEN
               ALLOCATE(scratch_sym(nnp%ang(ind)%symfgrp(s)%n_symf))
               scratch_max_sym = SIZE(scratch_sym)
            ELSE IF (scratch_max_sym < nnp%ang(ind)%symfgrp(s)%n_symf) THEN
               DEALLOCATE(scratch_sym)
               ALLOCATE(scratch_sym(nnp%ang(ind)%symfgrp(s)%n_symf))
               scratch_max_sym = SIZE(scratch_sym)
            END IF
            symtmp => scratch_sym(1:nnp%ang(ind)%symfgrp(s)%n_symf)
            !loop over associated neighbors
            IF (nnp%ang(ind)%symfgrp(s)%ele(1) == nnp%ang(ind)%symfgrp(s)%ele(2)) THEN
               DO j = 1, neighbor%n_ang1(s)
                  rvect1 = neighbor%dist_ang1(1:3, j, s)
                  r1 = neighbor%dist_ang1(4, j, s)
                  DO k = j + 1, neighbor%n_ang1(s)
                     rvect2 = neighbor%dist_ang1(1:3, k, s)
                     r2 = neighbor%dist_ang1(4, k, s)
                     rvect3 = rvect2 - rvect1
                     r3 = NORM2(rvect3(:))
                     IF (r3 < nnp%ang(ind)%symfgrp(s)%cutoff) THEN
                        CALL nnp_calc_ang(nnp, ind, s, rvect1, rvect2, rvect3, r1, r2, r3, symtmp)
                        DO sf = 1, nnp%ang(ind)%symfgrp(s)%n_symf
                           m = off + nnp%ang(ind)%symfgrp(s)%symf(sf)
                           nnp%ang(ind)%y(m - off) = nnp%ang(ind)%y(m - off) + symtmp(sf)
                        END DO
                     END IF
                  END DO
               END DO
            ELSE
               DO j = 1, neighbor%n_ang1(s)
                  rvect1 = neighbor%dist_ang1(1:3, j, s)
                  r1 = neighbor%dist_ang1(4, j, s)
                  DO k = 1, neighbor%n_ang2(s)
                     rvect2 = neighbor%dist_ang2(1:3, k, s)
                     r2 = neighbor%dist_ang2(4, k, s)
                     rvect3 = rvect2 - rvect1
                     r3 = NORM2(rvect3(:))
                     IF (r3 < nnp%ang(ind)%symfgrp(s)%cutoff) THEN
                        CALL nnp_calc_ang(nnp, ind, s, rvect1, rvect2, rvect3, r1, r2, r3, symtmp)
                        !loop over associated sym fncts
                        DO sf = 1, nnp%ang(ind)%symfgrp(s)%n_symf
                           m = off + nnp%ang(ind)%symfgrp(s)%symf(sf)
                           nnp%ang(ind)%y(m - off) = nnp%ang(ind)%y(m - off) + symtmp(sf)
                        END DO
                     END IF
                  END DO
               END DO
            END IF
            NULLIFY(symtmp)
         END DO
         CALL timestop(handle_sf)
      END IF

      !check extrapolation
      CALL nnp_check_extrapolation(nnp, ind)

      IF (PRESENT(dsymdxyz)) THEN
         IF (PRESENT(stress)) THEN
            IF (PRESENT(active_atoms) .AND. PRESENT(n_active)) THEN
               CALL nnp_scale_acsf(nnp, ind, dsymdxyz, stress, active_atoms, n_active)
            ELSE
               CALL nnp_scale_acsf(nnp, ind, dsymdxyz, stress)
            END IF
         ELSE
            IF (PRESENT(active_atoms) .AND. PRESENT(n_active)) THEN
               CALL nnp_scale_acsf(nnp, ind, dsymdxyz, active_atoms=active_atoms, n_active=n_active)
            ELSE
               CALL nnp_scale_acsf(nnp, ind, dsymdxyz)
            END IF
         END IF
      ELSE
         CALL nnp_scale_acsf(nnp, ind)
      END IF

      ! Persistent buffer — no release here; freed at program end via SAVE finalisation.
      CALL timestop(handle)

   END SUBROUTINE nnp_calc_acsf

! **************************************************************************************************
!> \brief Check if the nnp is extrapolating
!> \param nnp ...
!> \param ind ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_check_extrapolation(nnp, ind)

      TYPE(nnp_type), INTENT(INOUT)                      :: nnp
      INTEGER, INTENT(IN)                                :: ind

      REAL(KIND=dp), PARAMETER                           :: threshold = 0.0001_dp

      INTEGER                                            :: j
      LOGICAL                                            :: extrapolate

      extrapolate = nnp%output_expol

      DO j = 1, nnp%n_rad(ind)
         IF (nnp%rad(ind)%y(j) - &
             nnp%rad(ind)%loc_max(j) > threshold) THEN
            extrapolate = .TRUE.
         ELSE IF (-nnp%rad(ind)%y(j) + &
                  nnp%rad(ind)%loc_min(j) > threshold) THEN
            extrapolate = .TRUE.
         END IF
      END DO
      DO j = 1, nnp%n_ang(ind)
         IF (nnp%ang(ind)%y(j) - &
             nnp%ang(ind)%loc_max(j) > threshold) THEN
            extrapolate = .TRUE.
         ELSE IF (-nnp%ang(ind)%y(j) + &
                  nnp%ang(ind)%loc_min(j) > threshold) THEN
            extrapolate = .TRUE.
         END IF
      END DO

      nnp%output_expol = extrapolate

   END SUBROUTINE nnp_check_extrapolation

! **************************************************************************************************
!> \brief Scale and center symetry functions (and gradients)
!> \param nnp ...
!> \param ind ...
!> \param dsymdxyz ...
!> \param stress ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_scale_acsf(nnp, ind, dsymdxyz, stress, active_atoms, n_active)

      TYPE(nnp_type), INTENT(INOUT)                      :: nnp
      INTEGER, INTENT(IN)                                :: ind
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(OUT), &
         OPTIONAL                                        :: dsymdxyz, stress
      INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN)        :: active_atoms
      INTEGER, OPTIONAL, INTENT(IN)                      :: n_active

      INTEGER                                            :: idx, j, k, off

      IF (nnp%center_acsf) THEN
         DO j = 1, nnp%n_rad(ind)
            nnp%arc(ind)%layer(1)%node(j) = &
               (nnp%rad(ind)%y(j) - nnp%rad(ind)%loc_av(j))
         END DO
         off = nnp%n_rad(ind)
         DO j = 1, nnp%n_ang(ind)
            nnp%arc(ind)%layer(1)%node(j + off) = &
               (nnp%ang(ind)%y(j) - nnp%ang(ind)%loc_av(j))
         END DO

         IF (nnp%scale_acsf) THEN
            DO j = 1, nnp%n_rad(ind)
               nnp%arc(ind)%layer(1)%node(j) = &
                  nnp%arc(ind)%layer(1)%node(j)/ &
                  (nnp%rad(ind)%loc_max(j) - nnp%rad(ind)%loc_min(j))* &
                  (nnp%scmax - nnp%scmin) + nnp%scmin
            END DO
            off = nnp%n_rad(ind)
            DO j = 1, nnp%n_ang(ind)
               nnp%arc(ind)%layer(1)%node(j + off) = &
                  nnp%arc(ind)%layer(1)%node(j + off)/ &
                  (nnp%ang(ind)%loc_max(j) - nnp%ang(ind)%loc_min(j))* &
                  (nnp%scmax - nnp%scmin) + nnp%scmin
            END DO
         END IF
      ELSE IF (nnp%scale_acsf) THEN
         DO j = 1, nnp%n_rad(ind)
            nnp%arc(ind)%layer(1)%node(j) = &
               (nnp%rad(ind)%y(j) - nnp%rad(ind)%loc_min(j))/ &
               (nnp%rad(ind)%loc_max(j) - nnp%rad(ind)%loc_min(j))* &
               (nnp%scmax - nnp%scmin) + nnp%scmin
         END DO
         off = nnp%n_rad(ind)
         DO j = 1, nnp%n_ang(ind)
            nnp%arc(ind)%layer(1)%node(j + off) = &
               (nnp%ang(ind)%y(j) - nnp%ang(ind)%loc_min(j))/ &
               (nnp%ang(ind)%loc_max(j) - nnp%ang(ind)%loc_min(j))* &
               (nnp%scmax - nnp%scmin) + nnp%scmin
         END DO
      ELSE IF (nnp%scale_sigma_acsf) THEN
         DO j = 1, nnp%n_rad(ind)
            nnp%arc(ind)%layer(1)%node(j) = &
               (nnp%rad(ind)%y(j) - nnp%rad(ind)%loc_av(j))/ &
               nnp%rad(ind)%sigma(j)* &
               (nnp%scmax - nnp%scmin) + nnp%scmin
         END DO
         off = nnp%n_rad(ind)
         DO j = 1, nnp%n_ang(ind)
            nnp%arc(ind)%layer(1)%node(j + off) = &
               (nnp%ang(ind)%y(j) - nnp%ang(ind)%loc_av(j))/ &
               nnp%ang(ind)%sigma(j)* &
               (nnp%scmax - nnp%scmin) + nnp%scmin
         END DO
      ELSE
         DO j = 1, nnp%n_rad(ind)
            nnp%arc(ind)%layer(1)%node(j) = nnp%rad(ind)%y(j)
         END DO
         off = nnp%n_rad(ind)
         DO j = 1, nnp%n_ang(ind)
            nnp%arc(ind)%layer(1)%node(j + off) = nnp%ang(ind)%y(j)
         END DO
      END IF

      ! dsymdxyz scaling: walk only the active atoms when the caller provided them
      ! (drops O(num_atoms) per centre to O(n_active) ≈ O(1) at fixed density).
      IF (PRESENT(dsymdxyz)) THEN
         IF (nnp%scale_acsf) THEN
            IF (PRESENT(active_atoms) .AND. PRESENT(n_active)) THEN
               DO idx = 1, n_active
                  k = active_atoms(idx)
                  DO j = 1, nnp%n_rad(ind)
                     dsymdxyz(:, j, k) = dsymdxyz(:, j, k)/ &
                                         (nnp%rad(ind)%loc_max(j) - &
                                          nnp%rad(ind)%loc_min(j))* &
                                         (nnp%scmax - nnp%scmin)
                  END DO
               END DO
               off = nnp%n_rad(ind)
               DO idx = 1, n_active
                  k = active_atoms(idx)
                  DO j = 1, nnp%n_ang(ind)
                     dsymdxyz(:, j + off, k) = dsymdxyz(:, j + off, k)/ &
                                               (nnp%ang(ind)%loc_max(j) - &
                                                nnp%ang(ind)%loc_min(j))* &
                                               (nnp%scmax - nnp%scmin)
                  END DO
               END DO
            ELSE
               DO k = 1, nnp%num_atoms
                  DO j = 1, nnp%n_rad(ind)
                     dsymdxyz(:, j, k) = dsymdxyz(:, j, k)/ &
                                         (nnp%rad(ind)%loc_max(j) - &
                                          nnp%rad(ind)%loc_min(j))* &
                                         (nnp%scmax - nnp%scmin)
                  END DO
               END DO
               off = nnp%n_rad(ind)
               DO k = 1, nnp%num_atoms
                  DO j = 1, nnp%n_ang(ind)
                     dsymdxyz(:, j + off, k) = dsymdxyz(:, j + off, k)/ &
                                               (nnp%ang(ind)%loc_max(j) - &
                                                nnp%ang(ind)%loc_min(j))* &
                                               (nnp%scmax - nnp%scmin)
                  END DO
               END DO
            END IF
         ELSE IF (nnp%scale_sigma_acsf) THEN
            IF (PRESENT(active_atoms) .AND. PRESENT(n_active)) THEN
               DO idx = 1, n_active
                  k = active_atoms(idx)
                  DO j = 1, nnp%n_rad(ind)
                     dsymdxyz(:, j, k) = dsymdxyz(:, j, k)/ &
                                         nnp%rad(ind)%sigma(j)* &
                                         (nnp%scmax - nnp%scmin)
                  END DO
               END DO
               off = nnp%n_rad(ind)
               DO idx = 1, n_active
                  k = active_atoms(idx)
                  DO j = 1, nnp%n_ang(ind)
                     dsymdxyz(:, j + off, k) = dsymdxyz(:, j + off, k)/ &
                                               nnp%ang(ind)%sigma(j)* &
                                               (nnp%scmax - nnp%scmin)
                  END DO
               END DO
            ELSE
               DO k = 1, nnp%num_atoms
                  DO j = 1, nnp%n_rad(ind)
                     dsymdxyz(:, j, k) = dsymdxyz(:, j, k)/ &
                                         nnp%rad(ind)%sigma(j)* &
                                         (nnp%scmax - nnp%scmin)
                  END DO
               END DO
               off = nnp%n_rad(ind)
               DO k = 1, nnp%num_atoms
                  DO j = 1, nnp%n_ang(ind)
                     dsymdxyz(:, j + off, k) = dsymdxyz(:, j + off, k)/ &
                                               nnp%ang(ind)%sigma(j)* &
                                               (nnp%scmax - nnp%scmin)
                  END DO
               END DO
            END IF
         END IF
      END IF

      IF (PRESENT(stress)) THEN
         IF (nnp%scale_acsf) THEN
            DO j = 1, nnp%n_rad(ind)
               stress(:, :, j) = stress(:, :, j)/ &
                                 (nnp%rad(ind)%loc_max(j) - &
                                  nnp%rad(ind)%loc_min(j))* &
                                 (nnp%scmax - nnp%scmin)
            END DO
            off = nnp%n_rad(ind)
            DO j = 1, nnp%n_ang(ind)
               stress(:, :, j + off) = stress(:, :, j + off)/ &
                                       (nnp%ang(ind)%loc_max(j) - &
                                        nnp%ang(ind)%loc_min(j))* &
                                       (nnp%scmax - nnp%scmin)
            END DO
         ELSE IF (nnp%scale_sigma_acsf) THEN
            DO j = 1, nnp%n_rad(ind)
               stress(:, :, j) = stress(:, :, j)/ &
                                 nnp%rad(ind)%sigma(j)* &
                                 (nnp%scmax - nnp%scmin)
            END DO
            off = nnp%n_rad(ind)
            DO j = 1, nnp%n_ang(ind)
               stress(:, :, j + off) = stress(:, :, j + off)/ &
                                       nnp%ang(ind)%sigma(j)* &
                                       (nnp%scmax - nnp%scmin)
            END DO
         END IF
      END IF

   END SUBROUTINE nnp_scale_acsf

! **************************************************************************************************
!> \brief Calculate radial symmetry function and gradient (optinal)
!>        for given displacment vecotr rvect of atom i and j
!> \param nnp ...
!> \param ind ...
!> \param s ...
!> \param rvect ...
!> \param r ...
!> \param sym ...
!> \param force ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_calc_rad(nnp, ind, s, rvect, r, sym, force)

      TYPE(nnp_type), INTENT(IN)                         :: nnp
      INTEGER, INTENT(IN)                                :: ind, s
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rvect
      REAL(KIND=dp), INTENT(IN)                          :: r
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: sym
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT), &
         OPTIONAL                                        :: force

      INTEGER                                            :: k, sf, idx
      REAL(KIND=dp)                                      :: dfcutdr, dsymdr, eta, fcut, funccut, rs, &
                                                            tmp, t, cos_val, dcos_dx, tanh_val, rinv
      REAL(KIND=dp), DIMENSION(3)                        :: drdx

      !init
      drdx = 0.0_dp
      fcut = 0.0_dp
      dfcutdr = 0.0_dp

      !Calculate cutoff function and partial derivative
      funccut = nnp%rad(ind)%symfgrp(s)%cutoff !cutoff
      SELECT CASE (nnp%cut_type)
      CASE (nnp_cut_cos)
         tmp = r * nnp%rad(ind)%symfgrp(s)%pi_over_cut
         
         idx = INT((tmp - spline_cache%min_cos) * spline_cache%inv_dx_cos) + 1
         IF (idx < 1) idx = 1
         IF (idx > spline_cache%n_points - 1) idx = spline_cache%n_points - 1
         t = tmp - (spline_cache%min_cos + REAL(idx - 1, dp) * spline_cache%dx_cos)
         
         cos_val = spline_cache%coef_cos(1, idx) + t * (spline_cache%coef_cos(2, idx) + &
                   t * (spline_cache%coef_cos(3, idx) + t * spline_cache%coef_cos(4, idx)))
         
         fcut = 0.5_dp*(cos_val + 1.0_dp)
         
         IF (PRESENT(force)) THEN
            dcos_dx = spline_cache%coef_cos(2, idx) + t * (2.0_dp * spline_cache%coef_cos(3, idx) + &
                      t * 3.0_dp * spline_cache%coef_cos(4, idx))
            dfcutdr = 0.5_dp*(dcos_dx)*(pi/funccut)
         END IF

      CASE (nnp_cut_tanh)
         tmp = 1.0_dp - r * nnp%rad(ind)%symfgrp(s)%inv_cut
         
         idx = INT((tmp - spline_cache%min_tanh) * spline_cache%inv_dx_tanh) + 1
         IF (idx < 1) idx = 1
         IF (idx > spline_cache%n_points - 1) idx = spline_cache%n_points - 1
         t = tmp - (spline_cache%min_tanh + REAL(idx - 1, dp) * spline_cache%dx_tanh)
         
         tanh_val = spline_cache%coef_tanh(1, idx) + t * (spline_cache%coef_tanh(2, idx) + &
                    t * (spline_cache%coef_tanh(3, idx) + t * spline_cache%coef_tanh(4, idx)))
         
         fcut = tanh_val**3
         IF (PRESENT(force)) THEN
            dfcutdr = (-3.0_dp/funccut)*(tanh_val**2 - tanh_val**4)
         END IF

      CASE DEFAULT
         CPABORT("NNP| Cutoff function unknown")
      END SELECT

      ! 1 division -> 3 multiplies (division is ~6x slower on x86)
      IF (PRESENT(force)) THEN
         rinv = 1.0_dp/r
         drdx(1) = rvect(1)*rinv
         drdx(2) = rvect(2)*rinv
         drdx(3) = rvect(3)*rinv
      END IF

      !loop over sym fncts of sym fnct group s
      DO sf = 1, nnp%rad(ind)%symfgrp(s)%n_symf
         k = nnp%rad(ind)%symfgrp(s)%symf(sf) !symf indice
         eta = nnp%rad(ind)%eta(k) !eta
         rs = nnp%rad(ind)%rs(k) !rshift

         tmp = -eta*(r - rs)**2

         ! Fast Spline Evaluation for EXP(tmp).  We use sym(sf) (the spline value, ~EXP(tmp))
         ! directly in the chain rule instead of evaluating the spline derivative — the
         ! Hermite cubic value has O(h^4) error vs O(h^3) for its formal derivative, so this
         ! is both faster (skip 4 mul + 2 add) and slightly more accurate.
         IF (tmp <= spline_cache%min_exp) THEN
            sym(sf) = 0.0_dp
            IF (PRESENT(force)) dsymdr = 0.0_dp
         ELSE IF (tmp >= 0.0_dp) THEN
            sym(sf) = 1.0_dp
            IF (PRESENT(force)) dsymdr = 0.0_dp
         ELSE
            idx = INT((tmp - spline_cache%min_exp)*spline_cache%inv_dx_exp) + 1
            IF (idx < 1) idx = 1
            IF (idx > spline_cache%n_points - 1) idx = spline_cache%n_points - 1
            t = tmp - (spline_cache%min_exp + REAL(idx - 1, dp)*spline_cache%dx_exp)

            sym(sf) = spline_cache%coef_exp(1, idx) + t*(spline_cache%coef_exp(2, idx) + &
                      t*(spline_cache%coef_exp(3, idx) + t*spline_cache%coef_exp(4, idx)))

            IF (PRESENT(force)) dsymdr = sym(sf)*(-2.0_dp*eta*(r - rs))
         END IF

         ! combine radial symmetry function and cutoff function
         IF (PRESENT(force)) THEN
            force(:, sf) = fcut*dsymdr*drdx(:) + sym(sf)*dfcutdr*drdx(:)
         END IF
         sym(sf) = sym(sf)*fcut
      END DO

   END SUBROUTINE nnp_calc_rad

! **************************************************************************************************
!> \brief Calculate angular symmetry function and gradient (optinal)
!>        for given displacment vectors rvect1,2,3 of atom i,j and k
!> \param nnp ...
!> \param ind ...
!> \param s ...
!> \param rvect1 ...
!> \param rvect2 ...
!> \param rvect3 ...
!> \param r1 ...
!> \param r2 ...
!> \param r3 ...
!> \param sym ...
!> \param force ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_calc_ang(nnp, ind, s, rvect1, rvect2, rvect3, r1, r2, r3, sym, force)

      TYPE(nnp_type), INTENT(IN)                         :: nnp
      INTEGER, INTENT(IN)                                :: ind, s
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rvect1, rvect2, rvect3
      REAL(KIND=dp), INTENT(IN)                          :: r1, r2, r3
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: sym
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(OUT), &
         OPTIONAL                                        :: force

      INTEGER                                            :: i, m, sf, idx
      REAL(KIND=dp) :: angular, costheta, dfcutdr1, dfcutdr2, dfcutdr3, dsymdr1, dsymdr2, dsymdr3, &
         eta, f, fcut, fcut1, fcut2, fcut3, ftot, g, g_sq_inv, lam, pref, prefzeta, rsqr1, rsqr2, &
         rsqr3, r_sum_sq, r1inv, r2inv, r3inv, symtmp, tmp, tmp1, tmp2, tmp3, tmpzeta, zeta, t
      REAL(KIND=dp) :: dcos_dx1, dcos_dx2, dcos_dx3
      REAL(KIND=dp) :: tanh_val1, tanh_val2, tanh_val3
      REAL(KIND=dp), DIMENSION(3) :: dangulardx1, dangulardx2, dangulardx3, dcosthetadx1, &
         dcosthetadx2, dcosthetadx3, dfdx1, dfdx2, dfdx3, dgdx1, dgdx2, dgdx3, dr1dx, dr2dx, dr3dx

      rsqr1 = r1**2
      rsqr2 = r2**2
      rsqr3 = r3**2

      !init
      dangulardx1 = 0.0_dp
      dangulardx2 = 0.0_dp
      dangulardx3 = 0.0_dp
      dr1dx = 0.0_dp
      dr2dx = 0.0_dp
      dr3dx = 0.0_dp
      ftot = 0.0_dp
      dfcutdr1 = 0.0_dp
      dfcutdr2 = 0.0_dp
      dfcutdr3 = 0.0_dp

      ! Calculate cos(theta)
      f = (rsqr3 - rsqr1 - rsqr2)
      g = -2.0_dp*r1*r2
      costheta = f/g

      ! Calculate cutoff function and partial derivatives
      fcut = nnp%ang(ind)%symfgrp(s)%cutoff !cutoff
      SELECT CASE (nnp%cut_type)
      CASE (nnp_cut_cos)
         tmp1 = r1 * nnp%ang(ind)%symfgrp(s)%pi_over_cut
         idx = INT((tmp1 - spline_cache%min_cos) * spline_cache%inv_dx_cos) + 1
         IF (idx < 1) idx = 1; IF (idx > spline_cache%n_points - 1) idx = spline_cache%n_points - 1
         t = tmp1 - (spline_cache%min_cos + REAL(idx - 1, dp) * spline_cache%dx_cos)
         fcut1 = 0.5_dp * (spline_cache%coef_cos(1, idx) + t * (spline_cache%coef_cos(2, idx) + &
                 t * (spline_cache%coef_cos(3, idx) + t * spline_cache%coef_cos(4, idx))) + 1.0_dp)
         IF (PRESENT(force)) dcos_dx1 = spline_cache%coef_cos(2, idx) + & 
         t * (2.0_dp * spline_cache%coef_cos(3, idx) + t * 3.0_dp * spline_cache%coef_cos(4, idx))

         tmp2 = r2 * nnp%ang(ind)%symfgrp(s)%pi_over_cut
         idx = INT((tmp2 - spline_cache%min_cos) * spline_cache%inv_dx_cos) + 1
         IF (idx < 1) idx = 1; IF (idx > spline_cache%n_points - 1) idx = spline_cache%n_points - 1
         t = tmp2 - (spline_cache%min_cos + REAL(idx - 1, dp) * spline_cache%dx_cos)
         fcut2 = 0.5_dp * (spline_cache%coef_cos(1, idx) + t * (spline_cache%coef_cos(2, idx) + &
                 t * (spline_cache%coef_cos(3, idx) + t * spline_cache%coef_cos(4, idx))) + 1.0_dp)
         IF (PRESENT(force)) dcos_dx2 = spline_cache%coef_cos(2, idx) + &
               t * (2.0_dp * spline_cache%coef_cos(3, idx) + t * 3.0_dp * spline_cache%coef_cos(4, idx))

         tmp3 = r3 * nnp%ang(ind)%symfgrp(s)%pi_over_cut
         idx = INT((tmp3 - spline_cache%min_cos) * spline_cache%inv_dx_cos) + 1
         IF (idx < 1) idx = 1; IF (idx > spline_cache%n_points - 1) idx = spline_cache%n_points - 1
         t = tmp3 - (spline_cache%min_cos + REAL(idx - 1, dp) * spline_cache%dx_cos)
         fcut3 = 0.5_dp * (spline_cache%coef_cos(1, idx) + t * (spline_cache%coef_cos(2, idx) + &
                 t * (spline_cache%coef_cos(3, idx) + t * spline_cache%coef_cos(4, idx))) + 1.0_dp)
         IF (PRESENT(force)) dcos_dx3 = spline_cache%coef_cos(2, idx) + &
                  t * (2.0_dp * spline_cache%coef_cos(3, idx) + t * 3.0_dp * spline_cache%coef_cos(4, idx))

         ftot = fcut1*fcut2*fcut3
         IF (PRESENT(force)) THEN
            pref = 0.5_dp*(pi/fcut)
            dfcutdr1 = pref*(dcos_dx1)*fcut2*fcut3
            dfcutdr2 = pref*(dcos_dx2)*fcut1*fcut3
            dfcutdr3 = pref*(dcos_dx3)*fcut1*fcut2
         END IF

      CASE (nnp_cut_tanh)
         tmp1 = 1.0_dp - r1 * nnp%ang(ind)%symfgrp(s)%inv_cut
         idx = INT((tmp1 - spline_cache%min_tanh) * spline_cache%inv_dx_tanh) + 1
         IF (idx < 1) idx = 1; IF (idx > spline_cache%n_points - 1) idx = spline_cache%n_points - 1
         t = tmp1 - (spline_cache%min_tanh + REAL(idx - 1, dp) * spline_cache%dx_tanh)
         tanh_val1 = spline_cache%coef_tanh(1, idx) + t * (spline_cache%coef_tanh(2, idx) + &
                     t * (spline_cache%coef_tanh(3, idx) + t * spline_cache%coef_tanh(4, idx)))
         fcut1 = tanh_val1**3

         tmp2 = 1.0_dp - r2 * nnp%ang(ind)%symfgrp(s)%inv_cut
         idx = INT((tmp2 - spline_cache%min_tanh) * spline_cache%inv_dx_tanh) + 1
         IF (idx < 1) idx = 1; IF (idx > spline_cache%n_points - 1) idx = spline_cache%n_points - 1
         t = tmp2 - (spline_cache%min_tanh + REAL(idx - 1, dp) * spline_cache%dx_tanh)
         tanh_val2 = spline_cache%coef_tanh(1, idx) + t * (spline_cache%coef_tanh(2, idx) + &
                     t * (spline_cache%coef_tanh(3, idx) + t * spline_cache%coef_tanh(4, idx)))
         fcut2 = tanh_val2**3

         tmp3 = 1.0_dp - r3 * nnp%ang(ind)%symfgrp(s)%inv_cut
         idx = INT((tmp3 - spline_cache%min_tanh) * spline_cache%inv_dx_tanh) + 1
         IF (idx < 1) idx = 1; IF (idx > spline_cache%n_points - 1) idx = spline_cache%n_points - 1
         t = tmp3 - (spline_cache%min_tanh + REAL(idx - 1, dp) * spline_cache%dx_tanh)
         tanh_val3 = spline_cache%coef_tanh(1, idx) + t * (spline_cache%coef_tanh(2, idx) + &
                     t * (spline_cache%coef_tanh(3, idx) + t * spline_cache%coef_tanh(4, idx)))
         fcut3 = tanh_val3**3

         ftot = fcut1*fcut2*fcut3
         IF (PRESENT(force)) THEN
            pref = -3.0_dp/fcut
            dfcutdr1 = pref*(tanh_val1**2 - tanh_val1**4)*fcut2*fcut3
            dfcutdr2 = pref*(tanh_val2**2 - tanh_val2**4)*fcut1*fcut3
            dfcutdr3 = pref*(tanh_val3**2 - tanh_val3**4)*fcut1*fcut2
         END IF

      CASE DEFAULT
         CPABORT("NNP| Cutoff function unknown")
      END SELECT

      IF (PRESENT(force)) THEN
         ! 1 division -> 3 multiplies per pair-distance
         r1inv = 1.0_dp/r1
         r2inv = 1.0_dp/r2
         r3inv = 1.0_dp/r3
         dr1dx(1) = rvect1(1)*r1inv; dr1dx(2) = rvect1(2)*r1inv; dr1dx(3) = rvect1(3)*r1inv
         dr2dx(1) = rvect2(1)*r2inv; dr2dx(2) = rvect2(2)*r2inv; dr2dx(3) = rvect2(3)*r2inv
         dr3dx(1) = rvect3(1)*r3inv; dr3dx(2) = rvect3(2)*r3inv; dr3dx(3) = rvect3(3)*r3inv
      END IF

      ! Hoist invariants out of the sf loop:
      !  * r_sum_sq = r1^2 + r2^2 + r3^2 (the EXP argument is -eta*r_sum_sq, only eta varies per sf)
      !  * g_sq_inv = 1/g^2 (used in dangular/dr; g = -2 r1 r2 doesn't depend on sf)
      r_sum_sq = rsqr1 + rsqr2 + rsqr3
      g_sq_inv = 1.0_dp/(g*g)

      !loop over associated sym fncts
      DO sf = 1, nnp%ang(ind)%symfgrp(s)%n_symf
         m = nnp%ang(ind)%symfgrp(s)%symf(sf)
         lam = nnp%ang(ind)%lam(m) !lambda
         zeta = nnp%ang(ind)%zeta(m) !zeta
         prefzeta = nnp%ang(ind)%prefzeta(m) ! 2**(1-zeta)
         eta = nnp%ang(ind)%eta(m) !eta

         tmp = (1.0_dp + lam*costheta)

         IF (tmp <= 0.0_dp) THEN
            sym(sf) = 0.0_dp
            IF (PRESENT(force)) force(:, :, sf) = 0.0_dp
            CYCLE
         END IF

         ! Calculate angular symmetry function.  zeta_is_int / zeta_int are precomputed once
         ! at init from the (post-sort) zeta values; no NINT or FP-equality test in the hot loop.
         IF (nnp%ang(ind)%zeta_is_int(m)) THEN
            tmpzeta = tmp**(nnp%ang(ind)%zeta_int(m) - 1)
         ELSE
            tmpzeta = tmp**(zeta - 1.0_dp)
         END IF
         angular = tmpzeta*tmp

         ! Spline Evaluation for EXP.  Use symtmp (~EXP value) for the chain rule below
         ! instead of evaluating the spline derivative separately — same answer, slightly
         ! more accurate (O(h^4) value vs O(h^3) formal derivative).
         tmp = -eta*r_sum_sq
         IF (tmp <= spline_cache%min_exp) THEN
            symtmp = 0.0_dp
         ELSE IF (tmp >= 0.0_dp) THEN
            symtmp = 1.0_dp
         ELSE
            idx = INT((tmp - spline_cache%min_exp)*spline_cache%inv_dx_exp) + 1
            IF (idx < 1) idx = 1
            IF (idx > spline_cache%n_points - 1) idx = spline_cache%n_points - 1
            t = tmp - (spline_cache%min_exp + REAL(idx - 1, dp)*spline_cache%dx_exp)

            symtmp = spline_cache%coef_exp(1, idx) + t*(spline_cache%coef_exp(2, idx) + &
                     t*(spline_cache%coef_exp(3, idx) + t*spline_cache%coef_exp(4, idx)))
         END IF

         IF (PRESENT(force)) THEN
            pref = zeta*tmpzeta*g_sq_inv
            DO i = 1, 3
               dfdx1(i) = -2.0_dp*lam*(rvect1(i) + rvect2(i))
               dfdx2(i) = 2.0_dp*lam*(rvect3(i) + rvect1(i))
               dfdx3(i) = 2.0_dp*lam*(rvect2(i) - rvect3(i))

               tmp1 = 2.0_dp*r2*dr1dx(i)
               tmp2 = 2.0_dp*r1*dr2dx(i)
               dgdx1(i) = -(tmp1 + tmp2)
               dgdx2(i) = tmp1
               dgdx3(i) = tmp2

               dcosthetadx1(i) = dfdx1(i)*g - lam*f*dgdx1(i)
               dcosthetadx2(i) = dfdx2(i)*g - lam*f*dgdx2(i)
               dcosthetadx3(i) = dfdx3(i)*g - lam*f*dgdx3(i)

               dangulardx1(i) = pref*dcosthetadx1(i)
               dangulardx2(i) = pref*dcosthetadx2(i)
               dangulardx3(i) = pref*dcosthetadx3(i)
            END DO

            pref = prefzeta
            tmp = -2.0_dp*symtmp*eta
            dsymdr1 = tmp*r1
            dsymdr2 = tmp*r2
            dsymdr3 = tmp*r3

            tmp = pref*symtmp*ftot
            tmp1 = pref*angular*(ftot*dsymdr1 + dfcutdr1*symtmp)
            tmp2 = pref*angular*(ftot*dsymdr2 + dfcutdr2*symtmp)
            tmp3 = pref*angular*(ftot*dsymdr3 + dfcutdr3*symtmp)
            DO i = 1, 3
               force(i, 1, sf) = tmp*dangulardx1(i) + tmp1*dr1dx(i) + tmp2*dr2dx(i)
               force(i, 2, sf) = tmp*dangulardx2(i) - tmp1*dr1dx(i) + tmp3*dr3dx(i)
               force(i, 3, sf) = tmp*dangulardx3(i) - tmp2*dr2dx(i) - tmp3*dr3dx(i)
            END DO
         END IF

         pref = prefzeta
         sym(sf) = pref*angular*symtmp*ftot
      END DO

   END SUBROUTINE nnp_calc_ang

! **************************************************************************************************
!> \brief Sort element array according to atomic number
!> \param ele ...
!> \param nuc_ele ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_sort_ele(ele, nuc_ele)
      CHARACTER(len=2), DIMENSION(:), INTENT(INOUT)      :: ele
      INTEGER, DIMENSION(:), INTENT(INOUT)               :: nuc_ele

      CHARACTER(len=2)                                   :: tmp_ele
      INTEGER                                            :: i, j, loc, minimum, tmp_nuc_ele

      ! Determine atomic number
      DO i = 1, SIZE(ele)
         CALL get_ptable_info(ele(i), number=nuc_ele(i))
      END DO

      ! Sort both arrays
      DO i = 1, SIZE(ele) - 1
         minimum = nuc_ele(i)
         loc = i
         DO j = i + 1, SIZE(ele)
            IF (nuc_ele(j) < minimum) THEN
               loc = j
               minimum = nuc_ele(j)
            END IF
         END DO
         tmp_nuc_ele = nuc_ele(i)
         nuc_ele(i) = nuc_ele(loc)
         nuc_ele(loc) = tmp_nuc_ele

         tmp_ele = ele(i)
         ele(i) = ele(loc)
         ele(loc) = tmp_ele

      END DO

   END SUBROUTINE nnp_sort_ele

! **************************************************************************************************
!> \brief Sort symmetry functions according to different criteria
!> \param nnp ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_sort_acsf(nnp)
      TYPE(nnp_type), INTENT(INOUT)                      :: nnp

      INTEGER                                            :: i, j, k, loc

      ! First sort is according to symmetry function type
      ! This is done manually, since data structures are separate
      ! Note: Bubble sort is OK here, since small sort + special
      DO i = 1, nnp%n_ele
         ! sort by cutoff
         ! rad
         DO j = 1, nnp%n_rad(i) - 1
            loc = j
            DO k = j + 1, nnp%n_rad(i)
               IF (nnp%rad(i)%funccut(loc) > nnp%rad(i)%funccut(k)) THEN
                  loc = k
               END IF
            END DO
            ! swap symfnct
            CALL nnp_swaprad(nnp%rad(i), j, loc)
         END DO

         ! sort by eta
         ! rad
         DO j = 1, nnp%n_rad(i) - 1
            loc = j
            DO k = j + 1, nnp%n_rad(i)
               IF (nnp%rad(i)%funccut(loc) == nnp%rad(i)%funccut(k) .AND. &
                   nnp%rad(i)%eta(loc) > nnp%rad(i)%eta(k)) THEN
                  loc = k
               END IF
            END DO
            ! swap symfnct
            CALL nnp_swaprad(nnp%rad(i), j, loc)
         END DO

         ! sort by rshift
         ! rad
         DO j = 1, nnp%n_rad(i) - 1
            loc = j
            DO k = j + 1, nnp%n_rad(i)
               IF (nnp%rad(i)%funccut(loc) == nnp%rad(i)%funccut(k) .AND. &
                   nnp%rad(i)%eta(loc) == nnp%rad(i)%eta(k) .AND. &
                   nnp%rad(i)%rs(loc) > nnp%rad(i)%rs(k)) THEN
                  loc = k
               END IF
            END DO
            ! swap symfnct
            CALL nnp_swaprad(nnp%rad(i), j, loc)
         END DO

         ! sort by ele
         ! rad
         DO j = 1, nnp%n_rad(i) - 1
            loc = j
            DO k = j + 1, nnp%n_rad(i)
               IF (nnp%rad(i)%funccut(loc) == nnp%rad(i)%funccut(k) .AND. &
                   nnp%rad(i)%eta(loc) == nnp%rad(i)%eta(k) .AND. &
                   nnp%rad(i)%rs(loc) == nnp%rad(i)%rs(k) .AND. &
                   nnp%rad(i)%nuc_ele(loc) > nnp%rad(i)%nuc_ele(k)) THEN
                  loc = k
               END IF
            END DO
            ! swap symfnct
            CALL nnp_swaprad(nnp%rad(i), j, loc)
         END DO

         ! ang
         ! sort by cutoff
         DO j = 1, nnp%n_ang(i) - 1
            loc = j
            DO k = j + 1, nnp%n_ang(i)
               IF (nnp%ang(i)%funccut(loc) > nnp%ang(i)%funccut(k)) THEN
                  loc = k
               END IF
            END DO
            ! swap symfnct
            CALL nnp_swapang(nnp%ang(i), j, loc)
         END DO

         ! sort by eta
         ! ang
         DO j = 1, nnp%n_ang(i) - 1
            loc = j
            DO k = j + 1, nnp%n_ang(i)
               IF (nnp%ang(i)%funccut(loc) == nnp%ang(i)%funccut(k) .AND. &
                   nnp%ang(i)%eta(loc) > nnp%ang(i)%eta(k)) THEN
                  loc = k
               END IF
            END DO
            ! swap symfnct
            CALL nnp_swapang(nnp%ang(i), j, loc)
         END DO

         ! sort by zeta
         ! ang
         DO j = 1, nnp%n_ang(i) - 1
            loc = j
            DO k = j + 1, nnp%n_ang(i)
               IF (nnp%ang(i)%funccut(loc) == nnp%ang(i)%funccut(k) .AND. &
                   nnp%ang(i)%eta(loc) == nnp%ang(i)%eta(k) .AND. &
                   nnp%ang(i)%zeta(loc) > nnp%ang(i)%zeta(k)) THEN
                  loc = k
               END IF
            END DO
            ! swap symfnct
            CALL nnp_swapang(nnp%ang(i), j, loc)
         END DO

         ! sort by lambda
         ! ang
         DO j = 1, nnp%n_ang(i) - 1
            loc = j
            DO k = j + 1, nnp%n_ang(i)
               IF (nnp%ang(i)%funccut(loc) == nnp%ang(i)%funccut(k) .AND. &
                   nnp%ang(i)%eta(loc) == nnp%ang(i)%eta(k) .AND. &
                   nnp%ang(i)%zeta(loc) == nnp%ang(i)%zeta(k) .AND. &
                   nnp%ang(i)%lam(loc) > nnp%ang(i)%lam(k)) THEN
                  loc = k
               END IF
            END DO
            ! swap symfnct
            CALL nnp_swapang(nnp%ang(i), j, loc)
         END DO

         ! sort by ele
         ! ang, ele1
         DO j = 1, nnp%n_ang(i) - 1
            loc = j
            DO k = j + 1, nnp%n_ang(i)
               IF (nnp%ang(i)%funccut(loc) == nnp%ang(i)%funccut(k) .AND. &
                   nnp%ang(i)%eta(loc) == nnp%ang(i)%eta(k) .AND. &
                   nnp%ang(i)%zeta(loc) == nnp%ang(i)%zeta(k) .AND. &
                   nnp%ang(i)%lam(loc) == nnp%ang(i)%lam(k) .AND. &
                   nnp%ang(i)%nuc_ele1(loc) > nnp%ang(i)%nuc_ele1(k)) THEN
                  loc = k
               END IF
            END DO
            ! swap symfnct
            CALL nnp_swapang(nnp%ang(i), j, loc)
         END DO
         ! ang, ele2
         DO j = 1, nnp%n_ang(i) - 1
            loc = j
            DO k = j + 1, nnp%n_ang(i)
               IF (nnp%ang(i)%funccut(loc) == nnp%ang(i)%funccut(k) .AND. &
                   nnp%ang(i)%eta(loc) == nnp%ang(i)%eta(k) .AND. &
                   nnp%ang(i)%zeta(loc) == nnp%ang(i)%zeta(k) .AND. &
                   nnp%ang(i)%lam(loc) == nnp%ang(i)%lam(k) .AND. &
                   nnp%ang(i)%nuc_ele1(loc) == nnp%ang(i)%nuc_ele1(k) .AND. &
                   nnp%ang(i)%nuc_ele2(loc) > nnp%ang(i)%nuc_ele2(k)) THEN
                  loc = k
               END IF
            END DO
            ! swap symfnct
            CALL nnp_swapang(nnp%ang(i), j, loc)
         END DO
      END DO

   END SUBROUTINE nnp_sort_acsf

! **************************************************************************************************
!> \brief Swap two radial symmetry functions
!> \param rad ...
!> \param i ...
!> \param j ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_swaprad(rad, i, j)
      TYPE(nnp_acsf_rad_type), INTENT(INOUT)             :: rad
      INTEGER, INTENT(IN)                                :: i, j

      CHARACTER(len=2)                                   :: tmpc
      INTEGER                                            :: tmpi
      REAL(KIND=dp)                                      :: tmpr

      tmpr = rad%funccut(i)
      rad%funccut(i) = rad%funccut(j)
      rad%funccut(j) = tmpr

      tmpr = rad%eta(i)
      rad%eta(i) = rad%eta(j)
      rad%eta(j) = tmpr

      tmpr = rad%rs(i)
      rad%rs(i) = rad%rs(j)
      rad%rs(j) = tmpr

      tmpc = rad%ele(i)
      rad%ele(i) = rad%ele(j)
      rad%ele(j) = tmpc

      tmpi = rad%nuc_ele(i)
      rad%nuc_ele(i) = rad%nuc_ele(j)
      rad%nuc_ele(j) = tmpi

   END SUBROUTINE nnp_swaprad

! **************************************************************************************************
!> \brief Swap two angular symmetry functions
!> \param ang ...
!> \param i ...
!> \param j ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_swapang(ang, i, j)
      TYPE(nnp_acsf_ang_type), INTENT(INOUT)             :: ang
      INTEGER, INTENT(IN)                                :: i, j

      CHARACTER(len=2)                                   :: tmpc
      INTEGER                                            :: tmpi
      REAL(KIND=dp)                                      :: tmpr

      tmpr = ang%funccut(i)
      ang%funccut(i) = ang%funccut(j)
      ang%funccut(j) = tmpr

      tmpr = ang%eta(i)
      ang%eta(i) = ang%eta(j)
      ang%eta(j) = tmpr

      tmpr = ang%zeta(i)
      ang%zeta(i) = ang%zeta(j)
      ang%zeta(j) = tmpr

      tmpr = ang%prefzeta(i)
      ang%prefzeta(i) = ang%prefzeta(j)
      ang%prefzeta(j) = tmpr

      tmpr = ang%lam(i)
      ang%lam(i) = ang%lam(j)
      ang%lam(j) = tmpr

      tmpc = ang%ele1(i)
      ang%ele1(i) = ang%ele1(j)
      ang%ele1(j) = tmpc

      tmpi = ang%nuc_ele1(i)
      ang%nuc_ele1(i) = ang%nuc_ele1(j)
      ang%nuc_ele1(j) = tmpi

      tmpc = ang%ele2(i)
      ang%ele2(i) = ang%ele2(j)
      ang%ele2(j) = tmpc

      tmpi = ang%nuc_ele2(i)
      ang%nuc_ele2(i) = ang%nuc_ele2(j)
      ang%nuc_ele2(j) = tmpi

   END SUBROUTINE nnp_swapang


! **************************************************************************************************
!> \brief Initialize symmetry function groups and pre-calculate Hermite Splines
!> \param nnp ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_init_acsf_groups(nnp)

      TYPE(nnp_type), INTENT(INOUT)                      :: nnp

      INTEGER                                            :: ang, i, j, k, rad, s, zint
      INTEGER                                            :: col, pos
      REAL(KIND=dp)                                      :: funccut

      ! Variables for inline Spline Initialization
      REAL(KIND=dp)                                      :: x0, x1, y0, y1, m0, m1, h

      ! Precompute integer-zeta companions (NINT and is-integer flag) once.  Sort has already
      ! finalised the order of zeta(:) by the time this routine runs, so a single pass suffices.
      DO i = 1, nnp%n_ele
         DO s = 1, nnp%n_ang(i)
            zint = NINT(nnp%ang(i)%zeta(s))
            nnp%ang(i)%zeta_int(s) = zint
            nnp%ang(i)%zeta_is_int(s) = (REAL(zint, dp) == nnp%ang(i)%zeta(s))
         END DO
      END DO

      !find out how many symmetry functions groups are needed
      DO i = 1, nnp%n_ele
         nnp%rad(i)%n_symfgrp = 0
         nnp%ang(i)%n_symfgrp = 0
         !search radial symmetry functions
         DO j = 1, nnp%n_ele
            funccut = -1.0_dp
            DO s = 1, nnp%n_rad(i)
               IF (nnp%rad(i)%ele(s) == nnp%ele(j)) THEN
                  IF (ABS(nnp%rad(i)%funccut(s) - funccut) > 1.0e-5_dp) THEN
                     nnp%rad(i)%n_symfgrp = nnp%rad(i)%n_symfgrp + 1
                     funccut = nnp%rad(i)%funccut(s)
                  END IF
               END IF
            END DO
         END DO
         !search angular symmetry functions
         DO j = 1, nnp%n_ele
            DO k = j, nnp%n_ele
               funccut = -1.0_dp
               DO s = 1, nnp%n_ang(i)
                  IF ((nnp%ang(i)%ele1(s) == nnp%ele(j) .AND. &
                       nnp%ang(i)%ele2(s) == nnp%ele(k)) .OR. &
                      (nnp%ang(i)%ele1(s) == nnp%ele(k) .AND. &
                       nnp%ang(i)%ele2(s) == nnp%ele(j))) THEN
                     IF (ABS(nnp%ang(i)%funccut(s) - funccut) > 1.0e-5_dp) THEN
                        nnp%ang(i)%n_symfgrp = nnp%ang(i)%n_symfgrp + 1
                        funccut = nnp%ang(i)%funccut(s)
                     END IF
                  END IF
               END DO
            END DO
         END DO
      END DO

      !allocate memory for symmetry functions group
      DO i = 1, nnp%n_ele
         ALLOCATE (nnp%rad(i)%symfgrp(nnp%rad(i)%n_symfgrp))
         ALLOCATE (nnp%ang(i)%symfgrp(nnp%ang(i)%n_symfgrp))
         DO j = 1, nnp%rad(i)%n_symfgrp
            nnp%rad(i)%symfgrp(j)%n_symf = 0
         END DO
         DO j = 1, nnp%ang(i)%n_symfgrp
            nnp%ang(i)%symfgrp(j)%n_symf = 0
         END DO
      END DO

      !init symmetry functions group
      DO i = 1, nnp%n_ele
         rad = 0
         ang = 0
         DO j = 1, nnp%n_ele
            funccut = -1.0_dp
            DO s = 1, nnp%n_rad(i)
               IF (nnp%rad(i)%ele(s) == nnp%ele(j)) THEN
                  IF (ABS(nnp%rad(i)%funccut(s) - funccut) > 1.0e-5_dp) THEN
                     rad = rad + 1
                     funccut = nnp%rad(i)%funccut(s)
                     nnp%rad(i)%symfgrp(rad)%cutoff = funccut
                     nnp%rad(i)%symfgrp(rad)%cutoff_sq = funccut*funccut
                     nnp%rad(i)%symfgrp(rad)%inv_cut = 1.0_dp/funccut
                     nnp%rad(i)%symfgrp(rad)%pi_over_cut = pi/funccut
                     ALLOCATE (nnp%rad(i)%symfgrp(rad)%ele(1))
                     ALLOCATE (nnp%rad(i)%symfgrp(rad)%ele_ind(1))
                     nnp%rad(i)%symfgrp(rad)%ele(1) = nnp%ele(j)
                     nnp%rad(i)%symfgrp(rad)%ele_ind(1) = j
                  END IF
                  nnp%rad(i)%symfgrp(rad)%n_symf = nnp%rad(i)%symfgrp(rad)%n_symf + 1
               END IF
            END DO
         END DO
         DO j = 1, nnp%n_ele
            DO k = j, nnp%n_ele
               funccut = -1.0_dp
               DO s = 1, nnp%n_ang(i)
                  IF ((nnp%ang(i)%ele1(s) == nnp%ele(j) .AND. &
                       nnp%ang(i)%ele2(s) == nnp%ele(k)) .OR. &
                      (nnp%ang(i)%ele1(s) == nnp%ele(k) .AND. &
                       nnp%ang(i)%ele2(s) == nnp%ele(j))) THEN
                     IF (ABS(nnp%ang(i)%funccut(s) - funccut) > 1.0e-5_dp) THEN
                        ang = ang + 1
                        funccut = nnp%ang(i)%funccut(s)
                        nnp%ang(i)%symfgrp(ang)%cutoff = funccut
                        nnp%ang(i)%symfgrp(ang)%cutoff_sq = funccut*funccut
                        nnp%ang(i)%symfgrp(ang)%inv_cut = 1.0_dp/funccut
                        nnp%ang(i)%symfgrp(ang)%pi_over_cut = pi/funccut
                        ALLOCATE (nnp%ang(i)%symfgrp(ang)%ele(2))
                        ALLOCATE (nnp%ang(i)%symfgrp(ang)%ele_ind(2))
                        nnp%ang(i)%symfgrp(ang)%ele(1) = nnp%ele(j)
                        nnp%ang(i)%symfgrp(ang)%ele(2) = nnp%ele(k)
                        nnp%ang(i)%symfgrp(ang)%ele_ind(1) = j
                        nnp%ang(i)%symfgrp(ang)%ele_ind(2) = k
                     END IF
                     nnp%ang(i)%symfgrp(ang)%n_symf = nnp%ang(i)%symfgrp(ang)%n_symf + 1
                  END IF
               END DO
            END DO
         END DO
      END DO

      !add symmetry functions to associated groups
      DO i = 1, nnp%n_ele
         DO j = 1, nnp%rad(i)%n_symfgrp
            ALLOCATE (nnp%rad(i)%symfgrp(j)%symf(nnp%rad(i)%symfgrp(j)%n_symf))
            rad = 0
            DO s = 1, nnp%n_rad(i)
               IF (nnp%rad(i)%ele(s) == nnp%rad(i)%symfgrp(j)%ele(1)) THEN
                  IF (ABS(nnp%rad(i)%funccut(s) - nnp%rad(i)%symfgrp(j)%cutoff) <= 1.0e-5_dp) THEN
                     rad = rad + 1
                     nnp%rad(i)%symfgrp(j)%symf(rad) = s
                  END IF
               END IF
            END DO
         END DO
         DO j = 1, nnp%ang(i)%n_symfgrp
            ALLOCATE (nnp%ang(i)%symfgrp(j)%symf(nnp%ang(i)%symfgrp(j)%n_symf))
            ang = 0
            DO s = 1, nnp%n_ang(i)
               IF ((nnp%ang(i)%ele1(s) == nnp%ang(i)%symfgrp(j)%ele(1) .AND. &
                    nnp%ang(i)%ele2(s) == nnp%ang(i)%symfgrp(j)%ele(2)) .OR. &
                   (nnp%ang(i)%ele1(s) == nnp%ang(i)%symfgrp(j)%ele(2) .AND. &
                    nnp%ang(i)%ele2(s) == nnp%ang(i)%symfgrp(j)%ele(1))) THEN
                  IF (ABS(nnp%ang(i)%funccut(s) - nnp%ang(i)%symfgrp(j)%cutoff) <= 1.0e-5_dp) THEN
                     ang = ang + 1
                     nnp%ang(i)%symfgrp(j)%symf(ang) = s
                  END IF
               END IF
            END DO
         END DO

            ! Build element->group mapping for radial groups: ele_grp_count and ele_grp_list
            IF (nnp%rad(i)%n_symfgrp > 0) THEN
               ALLOCATE(nnp%rad(i)%ele_grp_count(nnp%n_ele))
               ALLOCATE(nnp%rad(i)%ele_grp_list(nnp%rad(i)%n_symfgrp, nnp%n_ele))
               nnp%rad(i)%ele_grp_count(:) = 0
               nnp%rad(i)%ele_grp_list(:,:) = 0
               DO j = 1, nnp%rad(i)%n_symfgrp
                  DO k = 1, SIZE(nnp%rad(i)%symfgrp(j)%ele_ind)
                     IF (nnp%rad(i)%symfgrp(j)%ele_ind(k) > 0) THEN
                        INTEGER :: col, pos
                        col = nnp%rad(i)%symfgrp(j)%ele_ind(k)
                        pos = nnp%rad(i)%ele_grp_count(col) + 1
                        nnp%rad(i)%ele_grp_list(pos, col) = j
                        nnp%rad(i)%ele_grp_count(col) = pos
                     END IF
                  END DO
               END DO
            END IF

            ! Build element->group mappings for angular groups: ele1 and ele2
            IF (nnp%ang(i)%n_symfgrp > 0) THEN
               ALLOCATE(nnp%ang(i)%ele1_grp_count(nnp%n_ele))
               ALLOCATE(nnp%ang(i)%ele1_grp_list(nnp%ang(i)%n_symfgrp, nnp%n_ele))
               ALLOCATE(nnp%ang(i)%ele2_grp_count(nnp%n_ele))
               ALLOCATE(nnp%ang(i)%ele2_grp_list(nnp%ang(i)%n_symfgrp, nnp%n_ele))
               nnp%ang(i)%ele1_grp_count(:) = 0
               nnp%ang(i)%ele2_grp_count(:) = 0
               nnp%ang(i)%ele1_grp_list(:,:) = 0
               nnp%ang(i)%ele2_grp_list(:,:) = 0
               DO j = 1, nnp%ang(i)%n_symfgrp
                  DO k = 1, SIZE(nnp%ang(i)%symfgrp(j)%ele_ind)
                     IF (nnp%ang(i)%symfgrp(j)%ele_ind(k) > 0) THEN
                        INTEGER :: col, pos
                        col = nnp%ang(i)%symfgrp(j)%ele_ind(k)
                        ! If this group lists ele_ind(1) or ele_ind(2) depends on construction
                        ! add group j to both ele1 and ele2 lists if matching
                        IF (nnp%ang(i)%symfgrp(j)%ele_ind(1) == col) THEN
                           pos = nnp%ang(i)%ele1_grp_count(col) + 1
                           nnp%ang(i)%ele1_grp_list(pos, col) = j
                           nnp%ang(i)%ele1_grp_count(col) = pos
                        END IF
                        IF (nnp%ang(i)%symfgrp(j)%ele_ind(2) == col) THEN
                           pos = nnp%ang(i)%ele2_grp_count(col) + 1
                           nnp%ang(i)%ele2_grp_list(pos, col) = j
                           nnp%ang(i)%ele2_grp_count(col) = pos
                        END IF
                     END IF
                  END DO
               END DO
            END IF
      END DO

      ! -------------------------------------------------------------------------
      ! Inline Hermite Spline Initialization
      ! -------------------------------------------------------------------------
      IF (.NOT. spline_cache%initialized) THEN
         ! Setup EXP spline
         spline_cache%min_exp = -50.0_dp
         spline_cache%dx_exp = ABS(spline_cache%min_exp) / REAL(spline_cache%n_points - 1, dp)
         spline_cache%inv_dx_exp = 1.0_dp / spline_cache%dx_exp
         ALLOCATE(spline_cache%coef_exp(4, spline_cache%n_points - 1))

         h = spline_cache%dx_exp
         DO i = 1, spline_cache%n_points - 1
            x0 = spline_cache%min_exp + REAL(i - 1, dp) * h
            x1 = x0 + h
            y0 = EXP(x0); m0 = EXP(x0)
            y1 = EXP(x1); m1 = EXP(x1)
            spline_cache%coef_exp(1, i) = y0
            spline_cache%coef_exp(2, i) = m0
            spline_cache%coef_exp(3, i) = (3.0_dp*(y1 - y0) / (h*h)) - ((2.0_dp*m0 + m1) / h)
            spline_cache%coef_exp(4, i) = (2.0_dp*(y0 - y1) / (h*h*h)) + ((m0 + m1) / (h*h))
         END DO

         ! Setup COS spline.  Cutoff usage is f_c(x) = 0.5*(cos(pi*x) + 1) for x = r/cutoff in
         ! [0, 1), so the spline argument tmp = pi*r/cutoff lives in [0, pi] only.  Using [0, pi]
         ! instead of [-pi, pi] doubles the resolution at the same n_points (~16x improvement
         ! in pointwise Hermite-cubic error, since error scales as h^4).
         spline_cache%min_cos = 0.0_dp
         spline_cache%max_cos = acos(-1.0_dp)
         spline_cache%dx_cos = (spline_cache%max_cos - spline_cache%min_cos) / REAL(spline_cache%n_points - 1, dp)
         spline_cache%inv_dx_cos = 1.0_dp / spline_cache%dx_cos
         ALLOCATE(spline_cache%coef_cos(4, spline_cache%n_points - 1))

         h = spline_cache%dx_cos
         DO i = 1, spline_cache%n_points - 1
            x0 = spline_cache%min_cos + REAL(i - 1, dp) * h
            x1 = x0 + h
            y0 = COS(x0); m0 = -SIN(x0)
            y1 = COS(x1); m1 = -SIN(x1)
            spline_cache%coef_cos(1, i) = y0
            spline_cache%coef_cos(2, i) = m0
            spline_cache%coef_cos(3, i) = (3.0_dp*(y1 - y0) / (h*h)) - ((2.0_dp*m0 + m1) / h)
            spline_cache%coef_cos(4, i) = (2.0_dp*(y0 - y1) / (h*h*h)) + ((m0 + m1) / (h*h))
         END DO

         ! Setup TANH spline.  Cutoff usage is f_c(x) = tanh^3(1 - x) for x = r/cutoff in [0, 1),
         ! so the spline argument tmp = 1 - r/cutoff lives in [0, 1] only.  Using [0, 1] instead
         ! of [-10, 10] gives ~20x finer resolution at the same n_points (~1.6e5 improvement in
         ! pointwise error, way past what the regtest tolerance demands).
         spline_cache%min_tanh = 0.0_dp
         spline_cache%max_tanh = 1.0_dp
         spline_cache%dx_tanh = (spline_cache%max_tanh - spline_cache%min_tanh) / REAL(spline_cache%n_points - 1, dp)
         spline_cache%inv_dx_tanh = 1.0_dp / spline_cache%dx_tanh
         ALLOCATE(spline_cache%coef_tanh(4, spline_cache%n_points - 1))

         h = spline_cache%dx_tanh
         DO i = 1, spline_cache%n_points - 1
            x0 = spline_cache%min_tanh + REAL(i - 1, dp) * h
            x1 = x0 + h
            y0 = TANH(x0); m0 = 1.0_dp - y0*y0
            y1 = TANH(x1); m1 = 1.0_dp - y1*y1
            spline_cache%coef_tanh(1, i) = y0
            spline_cache%coef_tanh(2, i) = m0
            spline_cache%coef_tanh(3, i) = (3.0_dp*(y1 - y0) / (h*h)) - ((2.0_dp*m0 + m1) / h)
            spline_cache%coef_tanh(4, i) = (2.0_dp*(y0 - y1) / (h*h*h)) + ((m0 + m1) / (h*h))
         END DO
         
         spline_cache%initialized = .TRUE.
      END IF

   END SUBROUTINE nnp_init_acsf_groups

! **************************************************************************************************
!> \brief Write symmetry function information
!> \param nnp ...
!> \param para_env ...
!> \param printtag ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_write_acsf(nnp, para_env, printtag)
      TYPE(nnp_type), INTENT(INOUT)                      :: nnp
      TYPE(mp_para_env_type), POINTER                    :: para_env
      CHARACTER(LEN=*), INTENT(IN)                       :: printtag

      CHARACTER(len=default_string_length)               :: my_label
      INTEGER                                            :: i, j, unit_nr
      TYPE(cp_logger_type), POINTER                      :: logger

      NULLIFY (logger)
      logger => cp_get_default_logger()

      my_label = TRIM(printtag)//"| "
      IF (para_env%is_source()) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger)
         WRITE (unit_nr, '(1X,A,1X,10(I2,1X))') TRIM(my_label)//" Activation functions:", nnp%actfnct(:)
         DO i = 1, nnp%n_ele
            WRITE (unit_nr, *) TRIM(my_label)//" short range atomic symmetry functions element "// &
               nnp%ele(i)//":"
            DO j = 1, nnp%n_rad(i)
               WRITE (unit_nr, '(1X,A,1X,I3,1X,A2,1X,I2,1X,A2,11X,3(F6.3,1X))') TRIM(my_label), j, nnp%ele(i), 2, &
                  nnp%rad(i)%ele(j), nnp%rad(i)%eta(j), &
                  nnp%rad(i)%rs(j), nnp%rad(i)%funccut(j)
            END DO
            DO j = 1, nnp%n_ang(i)
               WRITE (unit_nr, '(1X,A,1X,I3,1X,A2,1X,I2,2(1X,A2),1X,4(F6.3,1X))') &
                  TRIM(my_label), j, nnp%ele(i), 3, &
                  nnp%ang(i)%ele1(j), nnp%ang(i)%ele2(j), &
                  nnp%ang(i)%eta(j), nnp%ang(i)%lam(j), &
                  nnp%ang(i)%zeta(j), nnp%ang(i)%funccut(j)
            END DO
         END DO
      END IF

      RETURN

   END SUBROUTINE nnp_write_acsf

! **************************************************************************************************
!> \brief Build neighbour list for centre atom i by walking only the 27 bins around its current bin
!>        in the pre-expanded-image cell list.  O(k) per atom where k is the average per-bin image
!>        occupancy — independent of total system size.  No pbc() call in the inner loop: each image
!>        carries an explicit translation = hmat·shift, so the displacement is a direct subtraction.
!>        cell_list_cache must be populated by nnp_prepare_cell_list_cache before any call here.
!> \param nnp ...
!> \param neighbor ...
!> \param i ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_compute_neighbors(nnp, neighbor, i)

      TYPE(nnp_type), INTENT(INOUT), POINTER             :: nnp
      TYPE(nnp_neighbor_type), INTENT(INOUT)             :: neighbor
      INTEGER, INTENT(IN)                                :: i

      INTEGER                                            :: bx, by, bz, cell_idx, ind, img, j, &
                                 j_type, n_count, s, k
      INTEGER, DIMENSION(3)                              :: bin
      REAL(KIND=dp)                                      :: cx, cy, cz, drx, dry, drz, max_cut_sq, &
                                                            norm, r2

      ind = nnp%ele_ind(i)
      max_cut_sq = nnp%max_cut*nnp%max_cut

      ASSOCIATE (cl => cell_list_cache)
         ! Hoist centre coordinates as scalars (no array temporaries inside the bin walk)
         cx = cl%coord_primary(1, i)
         cy = cl%coord_primary(2, i)
         cz = cl%coord_primary(3, i)
         CALL nnp_cell_list_bin_from_position([cx, cy, cz], bin)

         DO bz = MAX(1, bin(3) - cl%bin_span(3)), MIN(cl%nbin(3), bin(3) + cl%bin_span(3))
         DO by = MAX(1, bin(2) - cl%bin_span(2)), MIN(cl%nbin(2), bin(2) + cl%bin_span(2))
         DO bx = MAX(1, bin(1) - cl%bin_span(1)), MIN(cl%nbin(1), bin(1) + cl%bin_span(1))
            cell_idx = nnp_cell_list_linear_index(bx, by, bz)
            img = cl%head(cell_idx)
            DO WHILE (img > 0)
               j = cl%image_atom(img)
               ! Skip self-image (j == i AND shift == 0); other periodic images of i are valid.
               IF (j == i .AND. cl%image_shift(1, img) == 0 .AND. &
                   cl%image_shift(2, img) == 0 .AND. cl%image_shift(3, img) == 0) THEN
                  img = cl%next(img)
                  CYCLE
               END IF
               ! Manual scalar dr (avoid Fortran array temporaries)
               drx = cx - (cl%coord_primary(1, j) + cl%image_translation(1, img))
               dry = cy - (cl%coord_primary(2, j) + cl%image_translation(2, img))
               drz = cz - (cl%coord_primary(3, j) + cl%image_translation(3, img))
               r2 = drx*drx + dry*dry + drz*drz
               ! Cheap reject before SQRT: anything past max_cut can't enter any group's list.
               IF (r2 >= max_cut_sq) THEN
                  img = cl%next(img)
                  CYCLE
               END IF
               ! Hoist neighbour element type once (saves repeated indirect loads in the cmp).
               j_type = nnp%ele_ind(j)
               ! lazy sqrt: compute norm only if at least one group accepts this image
               norm = -1.0_dp
               ! Radial groups: use per-element group mapping to iterate only relevant groups
               IF (ALLOCATED(nnp%rad(ind)%ele_grp_count)) THEN
                  n_count = 0
                  DO k = 1, nnp%rad(ind)%ele_grp_count(j_type)
                     s = nnp%rad(ind)%ele_grp_list(k, j_type)
                     IF (r2 < nnp%rad(ind)%symfgrp(s)%cutoff_sq) THEN
                        IF (norm < 0.0_dp) norm = SQRT(r2)
                        n_count = neighbor%n_rad(s) + 1
                        neighbor%n_rad(s) = n_count
                        neighbor%ind_rad(n_count, s) = j
                        neighbor%dist_rad(1, n_count, s) = drx
                        neighbor%dist_rad(2, n_count, s) = dry
                        neighbor%dist_rad(3, n_count, s) = drz
                        neighbor%dist_rad(4, n_count, s) = norm
                     END IF
                  END DO
               END IF
               ! Angular groups: check ele1 and ele2 mappings
               IF (ALLOCATED(nnp%ang(ind)%ele1_grp_count)) THEN
                  DO k = 1, nnp%ang(ind)%ele1_grp_count(j_type)
                     s = nnp%ang(ind)%ele1_grp_list(k, j_type)
                     IF (r2 < nnp%ang(ind)%symfgrp(s)%cutoff_sq) THEN
                        IF (norm < 0.0_dp) norm = SQRT(r2)
                        n_count = neighbor%n_ang1(s) + 1
                        neighbor%n_ang1(s) = n_count
                        neighbor%ind_ang1(n_count, s) = j
                        neighbor%dist_ang1(1, n_count, s) = drx
                        neighbor%dist_ang1(2, n_count, s) = dry
                        neighbor%dist_ang1(3, n_count, s) = drz
                        neighbor%dist_ang1(4, n_count, s) = norm
                     END IF
                  END DO
               END IF
               IF (ALLOCATED(nnp%ang(ind)%ele2_grp_count)) THEN
                  DO k = 1, nnp%ang(ind)%ele2_grp_count(j_type)
                     s = nnp%ang(ind)%ele2_grp_list(k, j_type)
                     IF (r2 < nnp%ang(ind)%symfgrp(s)%cutoff_sq) THEN
                        IF (norm < 0.0_dp) norm = SQRT(r2)
                        n_count = neighbor%n_ang2(s) + 1
                        neighbor%n_ang2(s) = n_count
                        neighbor%ind_ang2(n_count, s) = j
                        neighbor%dist_ang2(1, n_count, s) = drx
                        neighbor%dist_ang2(2, n_count, s) = dry
                        neighbor%dist_ang2(3, n_count, s) = drz
                        neighbor%dist_ang2(4, n_count, s) = norm
                     END IF
                  END DO
               END IF
               img = cl%next(img)
            END DO
         END DO
         END DO
         END DO
      END ASSOCIATE

   END SUBROUTINE nnp_compute_neighbors

! **************************************************************************************************
!> \brief Prepare (lazily rebuild) the module-level cell-list cache for this force evaluation.
!>        Must be called once per force eval (from nnp_calc_energy_force) before the per-atom
!>        nnp_compute_neighbors loop.  Cheap path is O(N): refresh PBC-wrapped coords and check the
!>        Verlet skin.  Full rebuild (O(N*(2p+1)^3)) is triggered only when an atom has moved more
!>        than 0.5*skin or when num_atoms / cell vectors / cutoff have changed.
!> \param nnp ...
!> \date   2026-04-24
!> \author Claudio Malvino
! **************************************************************************************************
   SUBROUTINE nnp_prepare_cell_list_cache(nnp)
      TYPE(nnp_type), INTENT(INOUT), POINTER             :: nnp

      CHARACTER(len=*), PARAMETER                        :: routineN = 'nnp_prepare_cell_list_cache'
      REAL(KIND=dp), PARAMETER                           :: default_skin = 1.0_dp, &
                                                            cell_tol = 1.0e-12_dp

      INTEGER                                            :: handle, j
      LOGICAL                                            :: need_rebuild
      REAL(KIND=dp)                                      :: dispsq, half_skin_sq, list_cutoff, &
                                                            max_disp_sq, skin
      REAL(KIND=dp), DIMENSION(3)                        :: dr
      TYPE(cell_type), POINTER                           :: cell

      CALL timeset(routineN, handle)

      NULLIFY (cell)
      CALL nnp_env_get(nnp_env=nnp, cell=cell)

      ! Ensure the persistent per-centre neighbour buffer is allocated to max capacity.
      ! Does nothing on subsequent steps if shape is unchanged; without this each centre
      ! atom would do 9 size-N ALLOCATEs in nnp_calc_acsf, giving O(N^2) malloc traffic.
      CALL nnp_neighbor_ensure_persistent(nnp)

      ! Skin: cap at 10% of the physical cutoff so the image set stays small.
      skin = MIN(default_skin, 0.1_dp*nnp%max_cut)
      list_cutoff = nnp%max_cut + skin

      ! Allocate primary-coord buffers when num_atoms changes (or first call)
      IF (.NOT. ALLOCATED(cell_list_cache%coord_primary) .OR. &
          cell_list_cache%num_atoms /= nnp%num_atoms) THEN
         IF (ALLOCATED(cell_list_cache%coord_primary)) DEALLOCATE (cell_list_cache%coord_primary)
         IF (ALLOCATED(cell_list_cache%ref_coord_primary)) DEALLOCATE (cell_list_cache%ref_coord_primary)
         ALLOCATE (cell_list_cache%coord_primary(3, nnp%num_atoms))
         ALLOCATE (cell_list_cache%ref_coord_primary(3, nnp%num_atoms))
         cell_list_cache%initialized = .FALSE.
      END IF

      ! Refresh primary atom positions every step (PBC-wrapped to the simulation cell)
      DO j = 1, nnp%num_atoms
         cell_list_cache%coord_primary(:, j) = pbc(nnp%coord(:, j), cell)
      END DO

      ! Decide whether a full rebuild is needed
      need_rebuild = .NOT. cell_list_cache%initialized
      IF (.NOT. need_rebuild) need_rebuild = (cell_list_cache%num_atoms /= nnp%num_atoms)
      IF (.NOT. need_rebuild) need_rebuild = (ABS(cell_list_cache%exact_cutoff - nnp%max_cut) > cell_tol)
      IF (.NOT. need_rebuild) need_rebuild = (ABS(cell_list_cache%verlet_skin - skin) > cell_tol)
      IF (.NOT. need_rebuild) need_rebuild = ANY(cell_list_cache%perd /= cell%perd)
      IF (.NOT. need_rebuild) need_rebuild = ANY(ABS(cell_list_cache%hmat - cell%hmat) > cell_tol)

      IF (.NOT. need_rebuild) THEN
         ! Verlet skin check: any atom moved more than 0.5*skin from its reference?
         half_skin_sq = (0.5_dp*skin)**2
         max_disp_sq = 0.0_dp
         DO j = 1, nnp%num_atoms
            dr(:) = cell_list_cache%coord_primary(:, j) - cell_list_cache%ref_coord_primary(:, j)
            ! minimum image difference — handles atoms wrapping across a box face
            dr = pbc(dr, cell)
            dispsq = dr(1)**2 + dr(2)**2 + dr(3)**2
            IF (dispsq > max_disp_sq) max_disp_sq = dispsq
            IF (max_disp_sq > half_skin_sq) THEN
               need_rebuild = .TRUE.
               EXIT
            END IF
         END DO
      END IF

      IF (need_rebuild) CALL nnp_build_cell_list_cache(nnp, cell, nnp%max_cut, list_cutoff, skin)

      CALL timestop(handle)

   END SUBROUTINE nnp_prepare_cell_list_cache

! **************************************************************************************************
!> \brief Rebuild the cell-list cache from scratch.  Snapshots reference coords, expands all
!>        (atom, integer PBC shift) image entries, computes the bounding box of all image
!>        positions (padded by half a skin), lays out a coarse grid with bin_width ≈ list_cutoff,
!>        and links each image into its bin via head/next.
!> \param nnp ...
!> \param cell ...
!> \param exact_cutoff ...
!> \param list_cutoff ...
!> \param skin ...
!> \date   2026-04-24
!> \author Claudio Malvino
! **************************************************************************************************
   SUBROUTINE nnp_build_cell_list_cache(nnp, cell, exact_cutoff, list_cutoff, skin)
      TYPE(nnp_type), INTENT(IN), POINTER                :: nnp
      TYPE(cell_type), INTENT(IN), POINTER               :: cell
      REAL(KIND=dp), INTENT(IN)                          :: exact_cutoff, list_cutoff, skin

      CHARACTER(len=*), PARAMETER                        :: routineN = 'nnp_build_cell_list_cache'
      REAL(KIND=dp), PARAMETER                           :: rebuild_eps = 1.0e-10_dp

      INTEGER                                            :: cell_idx, d, handle, img, j, &
                                                            n_cells_total, n_img, tx, ty, tz
      INTEGER, DIMENSION(3)                              :: bin
      REAL(KIND=dp), DIMENSION(3)                        :: ipos

      CALL timeset(routineN, handle)

      ASSOCIATE (cl => cell_list_cache)

         ! Snapshot rebuild-trigger state and reference coords
         cl%num_atoms = nnp%num_atoms
         cl%exact_cutoff = exact_cutoff
         cl%list_cutoff = list_cutoff
         cl%verlet_skin = skin
         cl%perd = cell%perd
         cl%hmat = cell%hmat
         cl%h_inv = cell%h_inv
         cl%orthorhombic = cell%orthorhombic
         cl%ref_coord_primary(:, :) = cl%coord_primary(:, :)

         ! Number of integer PBC shifts per dim that span list_cutoff (skin-inflated)
         CALL nnp_compute_pbc_copies(cl%image_copies, cell, list_cutoff)

         n_img = nnp%num_atoms* &
                 (2*cl%image_copies(1) + 1)* &
                 (2*cl%image_copies(2) + 1)* &
                 (2*cl%image_copies(3) + 1)
         cl%n_images = n_img

         IF (ALLOCATED(cl%image_atom)) DEALLOCATE (cl%image_atom)
         IF (ALLOCATED(cl%image_shift)) DEALLOCATE (cl%image_shift)
         IF (ALLOCATED(cl%image_translation)) DEALLOCATE (cl%image_translation)
         ALLOCATE (cl%image_atom(n_img))
         ALLOCATE (cl%image_shift(3, n_img))
         ALLOCATE (cl%image_translation(3, n_img))

         ! Build image table and bounding box of all image positions
         cl%lower(:) = HUGE(1.0_dp)
         cl%upper(:) = -HUGE(1.0_dp)
         img = 0
         DO j = 1, nnp%num_atoms
            DO tz = -cl%image_copies(3), cl%image_copies(3)
            DO ty = -cl%image_copies(2), cl%image_copies(2)
            DO tx = -cl%image_copies(1), cl%image_copies(1)
               img = img + 1
               cl%image_atom(img) = j
               cl%image_shift(1, img) = tx
               cl%image_shift(2, img) = ty
               cl%image_shift(3, img) = tz
               cl%image_translation(:, img) = cell%hmat(:, 1)*REAL(tx, dp) + &
                                              cell%hmat(:, 2)*REAL(ty, dp) + &
                                              cell%hmat(:, 3)*REAL(tz, dp)
               ipos(:) = cl%coord_primary(:, j) + cl%image_translation(:, img)
               DO d = 1, 3
                  IF (ipos(d) < cl%lower(d)) cl%lower(d) = ipos(d)
                  IF (ipos(d) > cl%upper(d)) cl%upper(d) = ipos(d)
               END DO
            END DO
            END DO
            END DO
         END DO

         ! Pad bounding box: half a skin on each side so atoms can drift before rebuild
         cl%lower(:) = cl%lower(:) - (0.5_dp*skin + rebuild_eps)
         cl%upper(:) = cl%upper(:) + (0.5_dp*skin + rebuild_eps)

         ! Coarse grid: bin_width ≈ list_cutoff, so search range is 1 (27-bin walk)
         DO d = 1, 3
            cl%nbin(d) = MAX(1, CEILING((cl%upper(d) - cl%lower(d))/list_cutoff))
            cl%bin_width(d) = (cl%upper(d) - cl%lower(d))/REAL(cl%nbin(d), dp)
            IF (cl%nbin(d) == 1) THEN
               cl%bin_span(d) = 0
            ELSE
               cl%bin_span(d) = MIN(cl%nbin(d) - 1, CEILING(list_cutoff/cl%bin_width(d)))
               IF (cl%bin_span(d) < 1) cl%bin_span(d) = 1
            END IF
         END DO

         n_cells_total = cl%nbin(1)*cl%nbin(2)*cl%nbin(3)
         cl%n_cells = n_cells_total

         IF (ALLOCATED(cl%head)) DEALLOCATE (cl%head)
         IF (ALLOCATED(cl%next)) DEALLOCATE (cl%next)
         ALLOCATE (cl%head(n_cells_total))
         ALLOCATE (cl%next(n_img))
         cl%head(:) = -1
         cl%next(:) = -1

         ! Link each image into its bin's chain (LIFO insertion at the head)
         DO img = 1, n_img
            j = cl%image_atom(img)
            ipos(:) = cl%coord_primary(:, j) + cl%image_translation(:, img)
            CALL nnp_cell_list_bin_from_position(ipos, bin)
            cell_idx = nnp_cell_list_linear_index(bin(1), bin(2), bin(3))
            cl%next(img) = cl%head(cell_idx)
            cl%head(cell_idx) = img
         END DO

         cl%initialized = .TRUE.

      END ASSOCIATE

      CALL timestop(handle)

   END SUBROUTINE nnp_build_cell_list_cache

! **************************************************************************************************
!> \brief Map a Cartesian position to its 1-based bin indices in the cell-list grid (clamped).
!> \param pos ...
!> \param bin ...
! **************************************************************************************************
   SUBROUTINE nnp_cell_list_bin_from_position(pos, bin)
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: pos
      INTEGER, DIMENSION(3), INTENT(OUT)                 :: bin

      INTEGER                                            :: d

      DO d = 1, 3
         bin(d) = INT((pos(d) - cell_list_cache%lower(d))/cell_list_cache%bin_width(d)) + 1
         IF (bin(d) < 1) bin(d) = 1
         IF (bin(d) > cell_list_cache%nbin(d)) bin(d) = cell_list_cache%nbin(d)
      END DO

   END SUBROUTINE nnp_cell_list_bin_from_position

! **************************************************************************************************
!> \brief Linearise 1-based 3D bin coordinates into the head() index.
!> \param bx ...
!> \param by ...
!> \param bz ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION nnp_cell_list_linear_index(bx, by, bz) RESULT(idx)
      INTEGER, INTENT(IN)                                :: bx, by, bz
      INTEGER                                            :: idx

      idx = bx + cell_list_cache%nbin(1)*((by - 1) + cell_list_cache%nbin(2)*(bz - 1))

   END FUNCTION nnp_cell_list_linear_index

! **************************************************************************************************
!> \brief Determine required pbc copies for small cells
!> \param pbc_copies ...
!> \param cell ...
!> \param cutoff ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_compute_pbc_copies(pbc_copies, cell, cutoff)
      INTEGER, DIMENSION(3), INTENT(INOUT)               :: pbc_copies
      TYPE(cell_type), INTENT(IN), POINTER               :: cell
      REAL(KIND=dp), INTENT(IN)                          :: cutoff

      REAL(KIND=dp)                                      :: proja, projb, projc
      REAL(KIND=dp), DIMENSION(3)                        :: axb, axc, bxc

      axb(1) = cell%hmat(2, 1)*cell%hmat(3, 2) - cell%hmat(3, 1)*cell%hmat(2, 2)
      axb(2) = cell%hmat(3, 1)*cell%hmat(1, 2) - cell%hmat(1, 1)*cell%hmat(3, 2)
      axb(3) = cell%hmat(1, 1)*cell%hmat(2, 2) - cell%hmat(2, 1)*cell%hmat(1, 2)
      axb(:) = axb(:)/NORM2(axb(:))

      axc(1) = cell%hmat(2, 1)*cell%hmat(3, 3) - cell%hmat(3, 1)*cell%hmat(2, 3)
      axc(2) = cell%hmat(3, 1)*cell%hmat(1, 3) - cell%hmat(1, 1)*cell%hmat(3, 3)
      axc(3) = cell%hmat(1, 1)*cell%hmat(2, 3) - cell%hmat(2, 1)*cell%hmat(1, 3)
      axc(:) = axc(:)/NORM2(axc(:))

      bxc(1) = cell%hmat(2, 2)*cell%hmat(3, 3) - cell%hmat(3, 2)*cell%hmat(2, 3)
      bxc(2) = cell%hmat(3, 2)*cell%hmat(1, 3) - cell%hmat(1, 2)*cell%hmat(3, 3)
      bxc(3) = cell%hmat(1, 2)*cell%hmat(2, 3) - cell%hmat(2, 2)*cell%hmat(1, 3)
      bxc(:) = bxc(:)/NORM2(bxc(:))

      proja = ABS(SUM(cell%hmat(:, 1)*bxc(:)))*0.5_dp
      projb = ABS(SUM(cell%hmat(:, 2)*axc(:)))*0.5_dp
      projc = ABS(SUM(cell%hmat(:, 3)*axb(:)))*0.5_dp

      pbc_copies(:) = 0
      DO WHILE ((pbc_copies(1) + 1)*proja <= cutoff)
         pbc_copies(1) = pbc_copies(1) + 1
      END DO
      DO WHILE ((pbc_copies(2) + 1)*projb <= cutoff)
         pbc_copies(2) = pbc_copies(2) + 1
      END DO
      DO WHILE ((pbc_copies(3) + 1)*projc <= cutoff)
         pbc_copies(3) = pbc_copies(3) + 1
      END DO
      ! Apply non periodic setting
      pbc_copies(:) = pbc_copies(:)*cell%perd(:)

   END SUBROUTINE nnp_compute_pbc_copies

! **************************************************************************************************
!> \brief Ensure the module-level persistent neighbour buffer is allocated to the maximum capacity
!>        needed across all element types.  After the first call (or after num_atoms / cutoff
!>        changes) every per-centre nnp_compute_neighbors invocation just resets counts and fills
!>        the existing buffer — no allocation traffic in the hot path.
!>        Capacity is sized for nnp%max_cut (consistent with the original buffer shape); the
!>        cell-list query never appends more neighbours than the per-group cutoff allows.
!> \param nnp ...
!> \date   2026-04-28
!> \author Claudio Malvino
! **************************************************************************************************
   SUBROUTINE nnp_neighbor_ensure_persistent(nnp)
      TYPE(nnp_type), INTENT(IN), POINTER                :: nnp

      CHARACTER(len=*), PARAMETER                        :: routineN = 'nnp_neighbor_ensure_persistent'

      INTEGER                                            :: handle, ie, n_ang_grp, n_capacity, &
                                 n_rad_grp, max_sym, s
      INTEGER, DIMENSION(3)                              :: pbc_copies
      TYPE(cell_type), POINTER                           :: cell

      CALL timeset(routineN, handle)

      NULLIFY (cell)
      CALL nnp_env_get(nnp_env=nnp, cell=cell)

      CALL nnp_compute_pbc_copies(pbc_copies, cell, nnp%max_cut)
      ! Correct capacity: product of copies per-dim (2*copies+1) times num_atoms
      n_capacity = nnp%num_atoms * (2*pbc_copies(1) + 1) * (2*pbc_copies(2) + 1) * (2*pbc_copies(3) + 1)

      n_rad_grp = 0
      n_ang_grp = 0
      DO ie = 1, nnp%n_ele
         IF (nnp%rad(ie)%n_symfgrp > n_rad_grp) n_rad_grp = nnp%rad(ie)%n_symfgrp
         IF (nnp%ang(ie)%n_symfgrp > n_ang_grp) n_ang_grp = nnp%ang(ie)%n_symfgrp
      END DO

      ! Ensure module-level scratch temporaries are allocated to handle the largest group
      ! Determine maximum symf size across all groups
      max_sym = 0
      DO ie = 1, nnp%n_ele
         DO s = 1, nnp%rad(ie)%n_symfgrp
            IF (nnp%rad(ie)%symfgrp(s)%n_symf > max_sym) max_sym = nnp%rad(ie)%symfgrp(s)%n_symf
         END DO
         DO s = 1, nnp%ang(ie)%n_symfgrp
            IF (nnp%ang(ie)%symfgrp(s)%n_symf > max_sym) max_sym = nnp%ang(ie)%symfgrp(s)%n_symf
         END DO
      END DO
      IF (max_sym < 1) max_sym = 1
      IF (.NOT. ALLOCATED(scratch_sym) .OR. scratch_max_sym < max_sym) THEN
         IF (ALLOCATED(scratch_sym)) DEALLOCATE(scratch_sym)
         IF (ALLOCATED(scratch_forcetmp)) DEALLOCATE(scratch_forcetmp)
         IF (ALLOCATED(scratch_force3tmp)) DEALLOCATE(scratch_force3tmp)
         ALLOCATE(scratch_sym(max_sym))
         ALLOCATE(scratch_forcetmp(3, max_sym))
         ALLOCATE(scratch_force3tmp(3,3,max_sym))
         scratch_max_sym = max_sym
      END IF

      ! Cheap path: buffer already at-or-above required shape, just record current pbc_copies
      IF (ALLOCATED(persistent_neighbor%dist_rad) .AND. &
          persistent_n_capacity >= n_capacity .AND. &
          persistent_n_rad_grp >= n_rad_grp .AND. &
          persistent_n_ang_grp >= n_ang_grp) THEN
         persistent_neighbor%pbc_copies = pbc_copies
         CALL timestop(handle)
         RETURN
      END IF

      ! (Re)allocate to max capacity across all element types
      IF (ALLOCATED(persistent_neighbor%dist_rad)) DEALLOCATE (persistent_neighbor%dist_rad)
      IF (ALLOCATED(persistent_neighbor%dist_ang1)) DEALLOCATE (persistent_neighbor%dist_ang1)
      IF (ALLOCATED(persistent_neighbor%dist_ang2)) DEALLOCATE (persistent_neighbor%dist_ang2)
      IF (ALLOCATED(persistent_neighbor%ind_rad)) DEALLOCATE (persistent_neighbor%ind_rad)
      IF (ALLOCATED(persistent_neighbor%ind_ang1)) DEALLOCATE (persistent_neighbor%ind_ang1)
      IF (ALLOCATED(persistent_neighbor%ind_ang2)) DEALLOCATE (persistent_neighbor%ind_ang2)
      IF (ALLOCATED(persistent_neighbor%n_rad)) DEALLOCATE (persistent_neighbor%n_rad)
      IF (ALLOCATED(persistent_neighbor%n_ang1)) DEALLOCATE (persistent_neighbor%n_ang1)
      IF (ALLOCATED(persistent_neighbor%n_ang2)) DEALLOCATE (persistent_neighbor%n_ang2)

      ALLOCATE (persistent_neighbor%dist_rad(4, n_capacity, n_rad_grp))
      ALLOCATE (persistent_neighbor%dist_ang1(4, n_capacity, n_ang_grp))
      ALLOCATE (persistent_neighbor%dist_ang2(4, n_capacity, n_ang_grp))
      ALLOCATE (persistent_neighbor%ind_rad(n_capacity, n_rad_grp))
      ALLOCATE (persistent_neighbor%ind_ang1(n_capacity, n_ang_grp))
      ALLOCATE (persistent_neighbor%ind_ang2(n_capacity, n_ang_grp))
      ALLOCATE (persistent_neighbor%n_rad(n_rad_grp))
      ALLOCATE (persistent_neighbor%n_ang1(n_ang_grp))
      ALLOCATE (persistent_neighbor%n_ang2(n_ang_grp))

      persistent_neighbor%pbc_copies = pbc_copies
      persistent_n_capacity = n_capacity
      persistent_n_rad_grp = n_rad_grp
      persistent_n_ang_grp = n_ang_grp

      CALL timestop(handle)

   END SUBROUTINE nnp_neighbor_ensure_persistent

END MODULE nnp_acsf
