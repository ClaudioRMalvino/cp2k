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
                                              nnp_arc_type,&
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
!> \brief Pre-expanded-image cell list used by the O(N) neighbour query.
!>        Each (atom j, integer PBC shift) is stored as an explicit image so
!>        the inner query loop does no pbc() calls.  Bin width is set so that
!>        each centre atom walks a 3x3x3 block of bins.  Rebuilds are
!>        triggered by Verlet-skin violation or by changes to num_atoms,
!>        cell vectors, or cutoff.
!> \date   2026-04-24
!> \author Claudio Malvino (crm98@cam.ac.uk)
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
      ! PBC-wrapped primary positions, refreshed every prepare call.
      REAL(KIND=dp), ALLOCATABLE           :: coord_primary(:, :)
      ! Snapshot of coord_primary at the last rebuild (Verlet reference).
      REAL(KIND=dp), ALLOCATABLE           :: ref_coord_primary(:, :)
      ! Image table: one entry per (atom j, integer PBC shift) pair.
      INTEGER, ALLOCATABLE                 :: image_atom(:)
      INTEGER, ALLOCATABLE                 :: image_shift(:, :)
      REAL(KIND=dp), ALLOCATABLE           :: image_translation(:, :)
      ! Linked-list cell grid.
      INTEGER, ALLOCATABLE                 :: head(:)
      INTEGER, ALLOCATABLE                 :: next(:)
      ! Longest bin chain at the last rebuild; bounds the per-centre
      ! candidate count (walk bins x occupancy) used to size the
      ! persistent neighbour buffer.
      INTEGER                              :: max_bin_occupancy = 0
   END TYPE nnp_cell_list_cache_type

   TYPE(nnp_cell_list_cache_type), SAVE, PRIVATE :: cell_list_cache

   ! Persistent per-centre neighbour buffer.  With centre-level OMP threading
   ! (nnp_calc_energy_force distributes whole centres over the team) every
   ! thread evaluates complete centres, so each thread owns its own buffer:
   ! all of the state below is THREADPRIVATE and grown lazily per thread in
   ! nnp_acsf_ensure_thread_buffers against the req_* sizes published once
   ! per force evaluation by nnp_neighbor_ensure_persistent (serial).
   TYPE(nnp_neighbor_type), SAVE, PRIVATE, TARGET :: persistent_neighbor
   INTEGER, SAVE, PRIVATE                          :: persistent_n_capacity = -1
   INTEGER, SAVE, PRIVATE                          :: persistent_n_rad_grp = -1
   INTEGER, SAVE, PRIVATE                          :: persistent_n_ang_grp = -1
   !$OMP THREADPRIVATE(persistent_neighbor, persistent_n_capacity, &
   !$OMP&              persistent_n_rad_grp, persistent_n_ang_grp)
   ! Required buffer sizes, written serially (outside any parallel region)
   ! once per force evaluation and read by all threads.
   INTEGER, SAVE, PRIVATE                          :: req_n_capacity = -1
   INTEGER, SAVE, PRIVATE                          :: req_n_rad_grp = -1
   INTEGER, SAVE, PRIVATE                          :: req_n_ang_grp = -1
   INTEGER, SAVE, PRIVATE                          :: req_max_sym = -1
   INTEGER, SAVE, PRIVATE                          :: req_max_rad = -1
   INTEGER, SAVE, PRIVATE                          :: req_max_ang = -1
   INTEGER, DIMENSION(3), SAVE, PRIVATE            :: req_pbc_copies = 0
   ! Per-thread scratch reused in the hot loops: per-pair outputs of
   ! nnp_calc_rad / nnp_calc_ang.  Accumulation into dsymdxyz / stress / y
   ! is direct (each thread owns the whole centre), so no atomics, critical
   ! sections, or per-group reduction buffers are needed.
   INTEGER, SAVE, PRIVATE                          :: scratch_max_sym = -1
   REAL(KIND=dp), ALLOCATABLE, SAVE, TARGET, PRIVATE :: scratch_sym(:)
   REAL(KIND=dp), ALLOCATABLE, SAVE, TARGET, PRIVATE :: scratch_forcetmp(:,:)
   REAL(KIND=dp), ALLOCATABLE, SAVE, TARGET, PRIVATE :: scratch_force3tmp(:,:,:)
   !$OMP THREADPRIVATE(scratch_max_sym, scratch_sym, scratch_forcetmp, scratch_force3tmp)
   ! Per-thread symmetry-function accumulators for the centre currently being
   ! evaluated.  These replace the shared nnp%rad(ind)%y / nnp%ang(ind)%y so
   ! two threads working on centres of the same element cannot collide.
   REAL(KIND=dp), ALLOCATABLE, SAVE, PRIVATE       :: ws_y_rad(:)
   REAL(KIND=dp), ALLOCATABLE, SAVE, PRIVATE       :: ws_y_ang(:)
   !$OMP THREADPRIVATE(ws_y_rad, ws_y_ang)
   ! Per-neighbour cutoff value/derivative and unit vector, all invariant
   ! over the inner k-loop of nnp_calc_ang.
   INTEGER, SAVE, PRIVATE                          :: scratch_max_fc = -1
   REAL(KIND=dp), ALLOCATABLE, SAVE, PRIVATE       :: scratch_fc1(:), scratch_dfc1(:)
   REAL(KIND=dp), ALLOCATABLE, SAVE, PRIVATE       :: scratch_fc2(:), scratch_dfc2(:)
   REAL(KIND=dp), ALLOCATABLE, SAVE, PRIVATE       :: scratch_dr1dx(:, :)  ! (3, n)
   REAL(KIND=dp), ALLOCATABLE, SAVE, PRIVATE       :: scratch_dr2dx(:, :)  ! (3, n)
   !$OMP THREADPRIVATE(scratch_max_fc, scratch_fc1, scratch_dfc1, scratch_fc2, scratch_dfc2, &
   !$OMP&              scratch_dr1dx, scratch_dr2dx)

   ! ------------------------------------------------------------------------------------------------
   ! Table-free minimax (Remez) polynomial replacements for the EXP/COS/TANH evaluations in the
   ! ACSF symmetry functions, evaluated by Clenshaw recurrence with compile-time-constant Chebyshev
   ! coefficients held in registers.  EXP uses range reduction (exp(x) = 2^k * exp(r), r in
   ! [-ln2/2, ln2/2]) so a degree-9 polynomial gives full RELATIVE precision over the whole
   ! physical range.  Validated against libm to < 1.4e-13 abs (cos/tanh) / 1.8e-14 rel (exp).
   ! The kernels live in THIS module (not splines_methods) so they share a translation unit with
   ! their only call sites: the FORCEINLINE directives are honoured without -ipo/-flto, the call
   ! overhead disappears, and the enclosing per-neighbour loops can vectorise.
   ! ------------------------------------------------------------------------------------------------
   REAL(KIND=dp), PARAMETER, PRIVATE :: ln2     = 0.69314718055994530942_dp
   REAL(KIND=dp), PARAMETER, PRIVATE :: inv_ln2 = 1.0_dp/0.69314718055994530942_dp

   ! cos on [0, pi]: degree 13 absolute minimax, Chebyshev basis
   REAL(KIND=dp), PARAMETER, PRIVATE, DIMENSION(0:13) :: c_cos = [ &
      2.88175278732944771e-15_dp, -1.13364817781174243e+00_dp, 5.57088022257035680e-15_dp, &
      1.38071776587197798e-01_dp, 5.42073505929545088e-15_dp, -4.49071424654963122e-03_dp, &
      5.43161196904855680e-15_dp, 6.77012758478224565e-05_dp, 5.75120890083155330e-15_dp, &
      -5.89129526975551519e-07_dp, 5.82145068847366206e-15_dp, 3.33806418687050215e-09_dp, &
      2.86158687413815466e-15_dp, -1.32954801664569218e-11_dp]

   ! tanh on [0, 1]: degree 14 absolute minimax, Chebyshev basis
   REAL(KIND=dp), PARAMETER, PRIVATE, DIMENSION(0:14) :: c_tanh = [ &
      4.20310759683689694e-01_dp, 3.83970922977466322e-01_dp, -4.06375122829675250e-02_dp, &
      -3.14432160508101854e-03_dp, 1.14609069487207785e-03_dp, -3.21132819936018998e-05_dp, &
      -2.25285844615948899e-05_dp, 2.67058666590674381e-06_dp, 2.67416778540631915e-07_dp, &
      -8.23968257221759303e-08_dp, 1.22045610095099853e-09_dp, 1.72044959173518527e-09_dp, &
      -1.76286593310703472e-10_dp, -2.28506334347878605e-11_dp, 5.88435334335686724e-12_dp]

   ! exp on reduced range [-ln2/2, ln2/2]: degree 9 relative minimax, Chebyshev basis
   REAL(KIND=dp), PARAMETER, PRIVATE, DIMENSION(0:9) :: c_exp = [ &
      1.03025449180961837e+00_dp, 3.51803207837704079e-01_dp, 3.03300103540964511e-02_dp, &
      1.74756361397688490e-03_dp, 7.55940398271240426e-05_dp, 2.61727190729881307e-06_dp, &
      7.55358006809067159e-08_dp, 1.86890629967635563e-09_dp, 4.04644375778992957e-11_dp, &
      7.74402774327636052e-13_dp]

CONTAINS

! **************************************************************************************************
!> \brief Calculate atom centered symmetry functions for given atom i
!>
!>        Thread-safe at centre level: when arc (and expol) are supplied, all
!>        mutable state touched by the call is either THREADPRIVATE module
!>        scratch or caller-owned (dsymdxyz, stress, active set, arc), so
!>        concurrent calls on different centres are safe.  Without arc/expol
!>        the shared nnp%arc(ind) / nnp%output_expol are used (serial callers
!>        such as the helium-NNP coupling keep their old behaviour).
!> \param nnp ...
!> \param i ...
!> \param dsymdxyz ...
!> \param stress ...
!> \param arc    (OPTIONAL) network state to fill; pass a per-thread copy when
!>               calling concurrently
!> \param expol  (OPTIONAL) extrapolation flag to update instead of
!>               nnp%output_expol
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_calc_acsf(nnp, i, dsymdxyz, stress, active_atoms, active_flag, n_active, &
                            arc, expol)

      TYPE(nnp_type), INTENT(INOUT), POINTER             :: nnp
      INTEGER, INTENT(IN)                                :: i
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT), &
         OPTIONAL                                        :: dsymdxyz
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT), &
         OPTIONAL, TARGET                                :: stress
      INTEGER, DIMENSION(:), OPTIONAL, INTENT(INOUT)     :: active_atoms, active_flag
      INTEGER, OPTIONAL, INTENT(OUT)                     :: n_active
      TYPE(nnp_arc_type), INTENT(INOUT), OPTIONAL, TARGET :: arc
      LOGICAL, INTENT(INOUT), OPTIONAL                   :: expol

      CHARACTER(len=*), PARAMETER                        :: routineN = 'nnp_calc_acsf'

      INTEGER                                            :: handle, handle_sf, ind, j, jj, k, kk, l, &
                                                            m, n_symf_s, off, s, sf, cut_type
      LOGICAL                                            :: has_stress, expol_loc
      REAL(KIND=dp)                                      :: r1, r2, r3, r3sq, r1inv, r2inv, &
                                                            fc_tmp, dfc_tmp, tanh_v, cos_v
      REAL(KIND=dp), DIMENSION(:), POINTER, CONTIGUOUS   :: symtmp
      REAL(KIND=dp), DIMENSION(:, :), POINTER, CONTIGUOUS :: forcetmp
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER, CONTIGUOUS :: force3tmp
      REAL(KIND=dp), DIMENSION(3)                        :: rvect1, rvect2, rvect3
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: stress_p
      TYPE(nnp_arc_type), POINTER                        :: arc_p
      TYPE(nnp_neighbor_type), POINTER                   :: neighbor

      CALL timeset(routineN, handle)

      !determine index of atom type
      ind = nnp%ele_ind(i)

      ! Cache cut_type locally; the compiler cannot prove the derived-type
      ! field is loop-invariant.
      cut_type = nnp%cut_type

      ! Grow this thread's persistent buffers to the sizes published by
      ! nnp_neighbor_ensure_persistent; a handful of integer compares on the
      ! fast path, allocations only on first use or system-size change.
      CALL nnp_acsf_ensure_thread_buffers()

      ! Reuse the thread's neighbour buffer; only the counts are reset
      ! between centres, so no allocations occur in the per-atom hot path.
      neighbor => persistent_neighbor
      neighbor%n_rad(:) = 0
      neighbor%n_ang1(:) = 0
      neighbor%n_ang2(:) = 0
      CALL nnp_compute_neighbors(nnp, neighbor, i)

      ! Network state to fill: the caller's per-thread copy if given, the
      ! shared instance otherwise.
      IF (PRESENT(arc)) THEN
         arc_p => arc
      ELSE
         arc_p => nnp%arc(ind)
      END IF

      has_stress = PRESENT(stress)
      IF (has_stress) THEN
         stress_p => stress
      ELSE
         NULLIFY(stress_p)
      END IF

      ! Build the active-atoms list (centre i plus every distinct neighbour
      ! referenced by the rad/ang lists).  active_flag is the dedup map; the
      ! caller resets touched entries to 0 between centres.
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
         ! Sort ascending so downstream walks hit dsymdxyz(:, :, k) in
         ! increasing k order; better prefetching and TLB locality.
         IF (n_active > 1) CALL nnp_sort_active_atoms(active_atoms, n_active)
      END IF

      ! Reset the per-thread symmetry-function accumulators:
      ws_y_rad(1:nnp%n_rad(ind)) = 0.0_dp
      ws_y_ang(1:nnp%n_ang(ind)) = 0.0_dp

      !calc forces
      IF (PRESENT(dsymdxyz)) THEN
         !loop over radial sym fnct grps
         CALL timeset('nnp_acsf_radial', handle_sf)
         DO s = 1, nnp%rad(ind)%n_symfgrp
            n_symf_s = nnp%rad(ind)%symfgrp(s)%n_symf
            symtmp => scratch_sym(1:n_symf_s)
            forcetmp => scratch_forcetmp(:, 1:n_symf_s)
            DO j = 1, neighbor%n_rad(s)
               rvect1 = neighbor%dist_rad(1:3, j, s)
               r1 = neighbor%dist_rad(4, j, s)
               CALL nnp_calc_rad(nnp, ind, s, rvect1, r1, symtmp, forcetmp)
               jj = neighbor%ind_rad(j, s)
               DO sf = 1, n_symf_s
                  m = nnp%rad(ind)%symfgrp(s)%symf(sf)
                  DO l = 1, 3
                     dsymdxyz(m, l, i) = dsymdxyz(m, l, i) + forcetmp(l, sf)
                     dsymdxyz(m, l, jj) = dsymdxyz(m, l, jj) - forcetmp(l, sf)
                  END DO
                  IF (has_stress) THEN
                     DO l = 1, 3
                        stress_p(1, l, m) = stress_p(1, l, m) + rvect1(1)*forcetmp(l, sf)
                        stress_p(2, l, m) = stress_p(2, l, m) + rvect1(2)*forcetmp(l, sf)
                        stress_p(3, l, m) = stress_p(3, l, m) + rvect1(3)*forcetmp(l, sf)
                     END DO
                  END IF
                  ws_y_rad(m) = ws_y_rad(m) + symtmp(sf)
               END DO
            END DO
         END DO
         NULLIFY(symtmp)
         NULLIFY(forcetmp)
         CALL timestop(handle_sf)
         !loop over angular sym fnct grps
         CALL timeset('nnp_acsf_angular', handle_sf)
         off = nnp%n_rad(ind)
         DO s = 1, nnp%ang(ind)%n_symfgrp
            n_symf_s = nnp%ang(ind)%symfgrp(s)%n_symf
            symtmp => scratch_sym(1:n_symf_s)
            force3tmp => scratch_force3tmp(:, :, 1:n_symf_s)

            ASSOCIATE (sg => nnp%ang(ind)%symfgrp(s))
            ! Per-j cutoff value, derivative, and unit vector are k-loop
            ! invariants; precompute once per group.  With the minimax
            ! kernels inlined (same translation unit) these loops vectorise.
            SELECT CASE (cut_type)
            CASE (nnp_cut_cos)
               DO j = 1, neighbor%n_ang1(s)
                  r1 = neighbor%dist_ang1(4, j, s)
                  fc_tmp = r1*sg%pi_over_cut
                  CALL minimax_cos_d(fc_tmp, cos_v, dfc_tmp)
                  scratch_fc1(j) = 0.5_dp*(cos_v + 1.0_dp)
                  scratch_dfc1(j) = 0.5_dp*dfc_tmp*sg%pi_over_cut
                  r1inv = 1.0_dp/r1
                  scratch_dr1dx(1, j) = neighbor%dist_ang1(1, j, s)*r1inv
                  scratch_dr1dx(2, j) = neighbor%dist_ang1(2, j, s)*r1inv
                  scratch_dr1dx(3, j) = neighbor%dist_ang1(3, j, s)*r1inv
               END DO
               IF (sg%ele(1) /= sg%ele(2)) THEN
                  DO k = 1, neighbor%n_ang2(s)
                     r2 = neighbor%dist_ang2(4, k, s)
                     fc_tmp = r2*sg%pi_over_cut
                     CALL minimax_cos_d(fc_tmp, cos_v, dfc_tmp)
                     scratch_fc2(k) = 0.5_dp*(cos_v + 1.0_dp)
                     scratch_dfc2(k) = 0.5_dp*dfc_tmp*sg%pi_over_cut
                     r2inv = 1.0_dp/r2
                     scratch_dr2dx(1, k) = neighbor%dist_ang2(1, k, s)*r2inv
                     scratch_dr2dx(2, k) = neighbor%dist_ang2(2, k, s)*r2inv
                     scratch_dr2dx(3, k) = neighbor%dist_ang2(3, k, s)*r2inv
                  END DO
               END IF
            CASE (nnp_cut_tanh)
               DO j = 1, neighbor%n_ang1(s)
                  r1 = neighbor%dist_ang1(4, j, s)
                  fc_tmp = 1.0_dp - r1*sg%inv_cut
                  tanh_v = minimax_tanh(fc_tmp)
                  scratch_fc1(j) = tanh_v**3
                  scratch_dfc1(j) = (-3.0_dp*sg%inv_cut)*(tanh_v**2 - tanh_v**4)
                  r1inv = 1.0_dp/r1
                  scratch_dr1dx(1, j) = neighbor%dist_ang1(1, j, s)*r1inv
                  scratch_dr1dx(2, j) = neighbor%dist_ang1(2, j, s)*r1inv
                  scratch_dr1dx(3, j) = neighbor%dist_ang1(3, j, s)*r1inv
               END DO
               IF (sg%ele(1) /= sg%ele(2)) THEN
                  DO k = 1, neighbor%n_ang2(s)
                     r2 = neighbor%dist_ang2(4, k, s)
                     fc_tmp = 1.0_dp - r2*sg%inv_cut
                     tanh_v = minimax_tanh(fc_tmp)
                     scratch_fc2(k) = tanh_v**3
                     scratch_dfc2(k) = (-3.0_dp*sg%inv_cut)*(tanh_v**2 - tanh_v**4)
                     r2inv = 1.0_dp/r2
                     scratch_dr2dx(1, k) = neighbor%dist_ang2(1, k, s)*r2inv
                     scratch_dr2dx(2, k) = neighbor%dist_ang2(2, k, s)*r2inv
                     scratch_dr2dx(3, k) = neighbor%dist_ang2(3, k, s)*r2inv
                  END DO
               END IF
            END SELECT

            IF (sg%ele(1) == sg%ele(2)) THEN
               DO j = 1, neighbor%n_ang1(s)
                  rvect1 = neighbor%dist_ang1(1:3, j, s)
                  r1 = neighbor%dist_ang1(4, j, s)
                  DO k = j + 1, neighbor%n_ang1(s)
                     rvect2 = neighbor%dist_ang1(1:3, k, s)
                     r2 = neighbor%dist_ang1(4, k, s)
                     rvect3 = rvect2 - rvect1
                     r3sq = rvect3(1)*rvect3(1) + rvect3(2)*rvect3(2) + rvect3(3)*rvect3(3)
                     IF (r3sq < sg%cutoff_sq) THEN
                        r3 = SQRT(r3sq)
                        jj = neighbor%ind_ang1(j, s)
                        kk = neighbor%ind_ang1(k, s)
                        CALL nnp_calc_ang(nnp, ind, s, rvect1, rvect2, rvect3, &
                                          r1, r2, r3, symtmp, force3tmp, &
                                          fc1_in=scratch_fc1(j), fc2_in=scratch_fc1(k), &
                                          dfc1dr_in=scratch_dfc1(j), dfc2dr_in=scratch_dfc1(k), &
                                          dr1dx_in=scratch_dr1dx(:, j), dr2dx_in=scratch_dr1dx(:, k))
                        DO sf = 1, n_symf_s
                           m = off + sg%symf(sf)
                           DO l = 1, 3
                              dsymdxyz(m, l, i) = dsymdxyz(m, l, i) + force3tmp(l, 1, sf)
                              dsymdxyz(m, l, jj) = dsymdxyz(m, l, jj) + force3tmp(l, 2, sf)
                              dsymdxyz(m, l, kk) = dsymdxyz(m, l, kk) + force3tmp(l, 3, sf)
                           END DO
                           IF (has_stress) THEN
                              DO l = 1, 3
                                 stress_p(1, l, m) = stress_p(1, l, m) &
                                    - rvect1(1)*force3tmp(l, 2, sf) - rvect2(1)*force3tmp(l, 3, sf)
                                 stress_p(2, l, m) = stress_p(2, l, m) &
                                    - rvect1(2)*force3tmp(l, 2, sf) - rvect2(2)*force3tmp(l, 3, sf)
                                 stress_p(3, l, m) = stress_p(3, l, m) &
                                    - rvect1(3)*force3tmp(l, 2, sf) - rvect2(3)*force3tmp(l, 3, sf)
                              END DO
                           END IF
                           ws_y_ang(sg%symf(sf)) = ws_y_ang(sg%symf(sf)) + symtmp(sf)
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
                     r3sq = rvect3(1)*rvect3(1) + rvect3(2)*rvect3(2) + rvect3(3)*rvect3(3)
                     IF (r3sq < sg%cutoff_sq) THEN
                        r3 = SQRT(r3sq)
                        jj = neighbor%ind_ang1(j, s)
                        kk = neighbor%ind_ang2(k, s)
                        CALL nnp_calc_ang(nnp, ind, s, rvect1, rvect2, rvect3, &
                                          r1, r2, r3, symtmp, force3tmp, &
                                          fc1_in=scratch_fc1(j), fc2_in=scratch_fc2(k), &
                                          dfc1dr_in=scratch_dfc1(j), dfc2dr_in=scratch_dfc2(k), &
                                          dr1dx_in=scratch_dr1dx(:, j), dr2dx_in=scratch_dr2dx(:, k))
                        DO sf = 1, n_symf_s
                           m = off + sg%symf(sf)
                           DO l = 1, 3
                              dsymdxyz(m, l, i) = dsymdxyz(m, l, i) + force3tmp(l, 1, sf)
                              dsymdxyz(m, l, jj) = dsymdxyz(m, l, jj) + force3tmp(l, 2, sf)
                              dsymdxyz(m, l, kk) = dsymdxyz(m, l, kk) + force3tmp(l, 3, sf)
                           END DO
                           IF (has_stress) THEN
                              DO l = 1, 3
                                 stress_p(1, l, m) = stress_p(1, l, m) &
                                    - rvect1(1)*force3tmp(l, 2, sf) - rvect2(1)*force3tmp(l, 3, sf)
                                 stress_p(2, l, m) = stress_p(2, l, m) &
                                    - rvect1(2)*force3tmp(l, 2, sf) - rvect2(2)*force3tmp(l, 3, sf)
                                 stress_p(3, l, m) = stress_p(3, l, m) &
                                    - rvect1(3)*force3tmp(l, 2, sf) - rvect2(3)*force3tmp(l, 3, sf)
                              END DO
                           END IF
                           ws_y_ang(sg%symf(sf)) = ws_y_ang(sg%symf(sf)) + symtmp(sf)
                        END DO
                     END IF
                  END DO
               END DO
            END IF
            END ASSOCIATE
         END DO
         NULLIFY(symtmp)
         NULLIFY(force3tmp)
         CALL timestop(handle_sf)
         !no forces
      ELSE
         !loop over radial sym fnct grps
         CALL timeset('nnp_acsf_radial', handle_sf)
         DO s = 1, nnp%rad(ind)%n_symfgrp
            n_symf_s = nnp%rad(ind)%symfgrp(s)%n_symf
            symtmp => scratch_sym(1:n_symf_s)
            DO j = 1, neighbor%n_rad(s)
               rvect1 = neighbor%dist_rad(1:3, j, s)
               r1 = neighbor%dist_rad(4, j, s)
               CALL nnp_calc_rad(nnp, ind, s, rvect1, r1, symtmp)
               DO sf = 1, n_symf_s
                  m = nnp%rad(ind)%symfgrp(s)%symf(sf)
                  ws_y_rad(m) = ws_y_rad(m) + symtmp(sf)
               END DO
            END DO
         END DO
         NULLIFY(symtmp)
         CALL timestop(handle_sf)
         !loop over angular sym fnct grps
         CALL timeset('nnp_acsf_angular', handle_sf)
         off = nnp%n_rad(ind)
         DO s = 1, nnp%ang(ind)%n_symfgrp
            n_symf_s = nnp%ang(ind)%symfgrp(s)%n_symf
            symtmp => scratch_sym(1:n_symf_s)

            ASSOCIATE (sg => nnp%ang(ind)%symfgrp(s))
            SELECT CASE (cut_type)
            CASE (nnp_cut_cos)
               DO j = 1, neighbor%n_ang1(s)
                  fc_tmp = neighbor%dist_ang1(4, j, s)*sg%pi_over_cut
                  scratch_fc1(j) = 0.5_dp*(minimax_cos(fc_tmp) + 1.0_dp)
               END DO
               IF (sg%ele(1) /= sg%ele(2)) THEN
                  DO k = 1, neighbor%n_ang2(s)
                     fc_tmp = neighbor%dist_ang2(4, k, s)*sg%pi_over_cut
                     scratch_fc2(k) = 0.5_dp*(minimax_cos(fc_tmp) + 1.0_dp)
                  END DO
               END IF
            CASE (nnp_cut_tanh)
               DO j = 1, neighbor%n_ang1(s)
                  fc_tmp = 1.0_dp - neighbor%dist_ang1(4, j, s)*sg%inv_cut
                  scratch_fc1(j) = minimax_tanh(fc_tmp)**3
               END DO
               IF (sg%ele(1) /= sg%ele(2)) THEN
                  DO k = 1, neighbor%n_ang2(s)
                     fc_tmp = 1.0_dp - neighbor%dist_ang2(4, k, s)*sg%inv_cut
                     scratch_fc2(k) = minimax_tanh(fc_tmp)**3
                  END DO
               END IF
            END SELECT

            IF (sg%ele(1) == sg%ele(2)) THEN
               DO j = 1, neighbor%n_ang1(s)
                  rvect1 = neighbor%dist_ang1(1:3, j, s)
                  r1 = neighbor%dist_ang1(4, j, s)
                  DO k = j + 1, neighbor%n_ang1(s)
                     rvect2 = neighbor%dist_ang1(1:3, k, s)
                     r2 = neighbor%dist_ang1(4, k, s)
                     rvect3 = rvect2 - rvect1
                     r3sq = rvect3(1)*rvect3(1) + rvect3(2)*rvect3(2) + rvect3(3)*rvect3(3)
                     IF (r3sq < sg%cutoff_sq) THEN
                        r3 = SQRT(r3sq)
                        CALL nnp_calc_ang(nnp, ind, s, rvect1, rvect2, rvect3, r1, r2, r3, symtmp, &
                                          fc1_in=scratch_fc1(j), fc2_in=scratch_fc1(k))
                        DO sf = 1, n_symf_s
                           ws_y_ang(sg%symf(sf)) = ws_y_ang(sg%symf(sf)) + symtmp(sf)
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
                     r3sq = rvect3(1)*rvect3(1) + rvect3(2)*rvect3(2) + rvect3(3)*rvect3(3)
                     IF (r3sq < sg%cutoff_sq) THEN
                        r3 = SQRT(r3sq)
                        CALL nnp_calc_ang(nnp, ind, s, rvect1, rvect2, rvect3, r1, r2, r3, symtmp, &
                                          fc1_in=scratch_fc1(j), fc2_in=scratch_fc2(k))
                        DO sf = 1, n_symf_s
                           ws_y_ang(sg%symf(sf)) = ws_y_ang(sg%symf(sf)) + symtmp(sf)
                        END DO
                     END IF
                  END DO
               END DO
            END IF
            END ASSOCIATE
         END DO
         NULLIFY(symtmp)
         CALL timestop(handle_sf)
      END IF

      !check extrapolation; update the caller's flag when given (thread-safe
      !path), the shared nnp%output_expol otherwise.
      IF (PRESENT(expol)) THEN
         CALL nnp_check_extrapolation(nnp, ind, expol)
      ELSE
         expol_loc = nnp%output_expol
         CALL nnp_check_extrapolation(nnp, ind, expol_loc)
         nnp%output_expol = expol_loc
      END IF

      IF (PRESENT(dsymdxyz)) THEN
         IF (PRESENT(stress)) THEN
            IF (PRESENT(active_atoms) .AND. PRESENT(n_active)) THEN
               CALL nnp_scale_acsf(nnp, ind, arc_p%layer(1)%node, dsymdxyz, stress, &
                                   active_atoms, n_active)
            ELSE
               CALL nnp_scale_acsf(nnp, ind, arc_p%layer(1)%node, dsymdxyz, stress)
            END IF
         ELSE
            IF (PRESENT(active_atoms) .AND. PRESENT(n_active)) THEN
               CALL nnp_scale_acsf(nnp, ind, arc_p%layer(1)%node, dsymdxyz, &
                                   active_atoms=active_atoms, n_active=n_active)
            ELSE
               CALL nnp_scale_acsf(nnp, ind, arc_p%layer(1)%node, dsymdxyz)
            END IF
         END IF
      ELSE
         CALL nnp_scale_acsf(nnp, ind, arc_p%layer(1)%node)
      END IF

      ! Module-level buffers are SAVE; no per-call release.
      CALL timestop(handle)

   END SUBROUTINE nnp_calc_acsf

! **************************************************************************************************
!> \brief Check if the nnp is extrapolating (reads the calling thread's
!>        ws_y_rad / ws_y_ang accumulators)
!> \param nnp ...
!> \param ind ...
!> \param expol  flag to OR the extrapolation status into
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_check_extrapolation(nnp, ind, expol)

      TYPE(nnp_type), INTENT(IN)                         :: nnp
      INTEGER, INTENT(IN)                                :: ind
      LOGICAL, INTENT(INOUT)                             :: expol

      REAL(KIND=dp), PARAMETER                           :: threshold = 0.0001_dp

      INTEGER                                            :: j
      LOGICAL                                            :: extrapolate

      extrapolate = expol

      DO j = 1, nnp%n_rad(ind)
         IF (ws_y_rad(j) - &
             nnp%rad(ind)%loc_max(j) > threshold) THEN
            extrapolate = .TRUE.
         ELSE IF (-ws_y_rad(j) + &
                  nnp%rad(ind)%loc_min(j) > threshold) THEN
            extrapolate = .TRUE.
         END IF
      END DO
      DO j = 1, nnp%n_ang(ind)
         IF (ws_y_ang(j) - &
             nnp%ang(ind)%loc_max(j) > threshold) THEN
            extrapolate = .TRUE.
         ELSE IF (-ws_y_ang(j) + &
                  nnp%ang(ind)%loc_min(j) > threshold) THEN
            extrapolate = .TRUE.
         END IF
      END DO

      expol = extrapolate

   END SUBROUTINE nnp_check_extrapolation

! **************************************************************************************************
!> \brief Scale and center symetry functions (and gradients).  Reads the
!>        calling thread's ws_y_rad / ws_y_ang accumulators and writes the
!>        caller-supplied input layer.
!> \param nnp ...
!> \param ind ...
!> \param node1  input-layer node vector to fill (layer 1 of the network
!>               state the caller will evaluate)
!> \param dsymdxyz ...
!> \param stress ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_scale_acsf(nnp, ind, node1, dsymdxyz, stress, active_atoms, n_active)

      TYPE(nnp_type), INTENT(INOUT)                      :: nnp
      INTEGER, INTENT(IN)                                :: ind
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: node1
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT), &
         OPTIONAL                                        :: dsymdxyz, stress
      INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN)        :: active_atoms
      INTEGER, OPTIONAL, INTENT(IN)                      :: n_active

      INTEGER                                            :: idx, j, k, off

      ! Build the input layer from the centre's symmetry-function values.
      ! Reciprocal scale factors are precomputed in nnp_init_model, so the
      ! per-symmetry-function divides become multiplies.
      IF (nnp%center_acsf) THEN
         DO j = 1, nnp%n_rad(ind)
            node1(j) = (ws_y_rad(j) - nnp%rad(ind)%loc_av(j))
         END DO
         off = nnp%n_rad(ind)
         DO j = 1, nnp%n_ang(ind)
            node1(j + off) = (ws_y_ang(j) - nnp%ang(ind)%loc_av(j))
         END DO

         IF (nnp%scale_acsf) THEN
            DO j = 1, nnp%n_rad(ind)
               node1(j) = node1(j)*nnp%rad(ind)%scale_factor(j) + nnp%scmin
            END DO
            off = nnp%n_rad(ind)
            DO j = 1, nnp%n_ang(ind)
               node1(j + off) = node1(j + off)*nnp%ang(ind)%scale_factor(j) + nnp%scmin
            END DO
         END IF
      ELSE IF (nnp%scale_acsf) THEN
         DO j = 1, nnp%n_rad(ind)
            node1(j) = (ws_y_rad(j) - nnp%rad(ind)%loc_min(j))*nnp%rad(ind)%scale_factor(j) &
                       + nnp%scmin
         END DO
         off = nnp%n_rad(ind)
         DO j = 1, nnp%n_ang(ind)
            node1(j + off) = (ws_y_ang(j) - nnp%ang(ind)%loc_min(j))*nnp%ang(ind)%scale_factor(j) &
                             + nnp%scmin
         END DO
      ELSE IF (nnp%scale_sigma_acsf) THEN
         DO j = 1, nnp%n_rad(ind)
            node1(j) = (ws_y_rad(j) - nnp%rad(ind)%loc_av(j))*nnp%rad(ind)%inv_sigma_factor(j) &
                       + nnp%scmin
         END DO
         off = nnp%n_rad(ind)
         DO j = 1, nnp%n_ang(ind)
            node1(j + off) = (ws_y_ang(j) - nnp%ang(ind)%loc_av(j))*nnp%ang(ind)%inv_sigma_factor(j) &
                             + nnp%scmin
         END DO
      ELSE
         DO j = 1, nnp%n_rad(ind)
            node1(j) = ws_y_rad(j)
         END DO
         off = nnp%n_rad(ind)
         DO j = 1, nnp%n_ang(ind)
            node1(j + off) = ws_y_ang(j)
         END DO
      END IF

      ! Scale dsymdxyz across the active atoms (serial: with centre-level
      ! threading this whole routine runs once per centre per thread).
      IF (PRESENT(dsymdxyz)) THEN
         IF (nnp%scale_acsf) THEN
            IF (PRESENT(active_atoms) .AND. PRESENT(n_active)) THEN
               DO idx = 1, n_active
                  k = active_atoms(idx)
                  DO j = 1, nnp%n_rad(ind)
                     dsymdxyz(j, :, k) = dsymdxyz(j, :, k)*nnp%rad(ind)%scale_factor(j)
                  END DO
               END DO
               off = nnp%n_rad(ind)
               DO idx = 1, n_active
                  k = active_atoms(idx)
                  DO j = 1, nnp%n_ang(ind)
                     dsymdxyz(j + off, :, k) = dsymdxyz(j + off, :, k)*nnp%ang(ind)%scale_factor(j)
                  END DO
               END DO
            ELSE
               DO k = 1, nnp%num_atoms
                  DO j = 1, nnp%n_rad(ind)
                     dsymdxyz(j, :, k) = dsymdxyz(j, :, k)*nnp%rad(ind)%scale_factor(j)
                  END DO
               END DO
               off = nnp%n_rad(ind)
               DO k = 1, nnp%num_atoms
                  DO j = 1, nnp%n_ang(ind)
                     dsymdxyz(j + off, :, k) = dsymdxyz(j + off, :, k)*nnp%ang(ind)%scale_factor(j)
                  END DO
               END DO
            END IF
         ELSE IF (nnp%scale_sigma_acsf) THEN
            IF (PRESENT(active_atoms) .AND. PRESENT(n_active)) THEN
               DO idx = 1, n_active
                  k = active_atoms(idx)
                  DO j = 1, nnp%n_rad(ind)
                     dsymdxyz(j, :, k) = dsymdxyz(j, :, k)*nnp%rad(ind)%inv_sigma_factor(j)
                  END DO
               END DO
               off = nnp%n_rad(ind)
               DO idx = 1, n_active
                  k = active_atoms(idx)
                  DO j = 1, nnp%n_ang(ind)
                     dsymdxyz(j + off, :, k) = dsymdxyz(j + off, :, k)*nnp%ang(ind)%inv_sigma_factor(j)
                  END DO
               END DO
            ELSE
               DO k = 1, nnp%num_atoms
                  DO j = 1, nnp%n_rad(ind)
                     dsymdxyz(j, :, k) = dsymdxyz(j, :, k)*nnp%rad(ind)%inv_sigma_factor(j)
                  END DO
               END DO
               off = nnp%n_rad(ind)
               DO k = 1, nnp%num_atoms
                  DO j = 1, nnp%n_ang(ind)
                     dsymdxyz(j + off, :, k) = dsymdxyz(j + off, :, k)*nnp%ang(ind)%inv_sigma_factor(j)
                  END DO
               END DO
            END IF
         END IF
      END IF

      IF (PRESENT(stress)) THEN
         IF (nnp%scale_acsf) THEN
            DO j = 1, nnp%n_rad(ind)
               stress(:, :, j) = stress(:, :, j)*nnp%rad(ind)%scale_factor(j)
            END DO
            off = nnp%n_rad(ind)
            DO j = 1, nnp%n_ang(ind)
               stress(:, :, j + off) = stress(:, :, j + off)*nnp%ang(ind)%scale_factor(j)
            END DO
         ELSE IF (nnp%scale_sigma_acsf) THEN
            DO j = 1, nnp%n_rad(ind)
               stress(:, :, j) = stress(:, :, j)*nnp%rad(ind)%inv_sigma_factor(j)
            END DO
            off = nnp%n_rad(ind)
            DO j = 1, nnp%n_ang(ind)
               stress(:, :, j + off) = stress(:, :, j + off)*nnp%ang(ind)%inv_sigma_factor(j)
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

      INTEGER                                            :: k, sf, cut_type
      REAL(KIND=dp)                                      :: dfcutdr, dsymdr, eta, fcut, funccut, rs, &
                                                            tmp, cos_val, dcos_dx, tanh_val, rinv
      REAL(KIND=dp), DIMENSION(3)                        :: drdx

      !init
      drdx = 0.0_dp
      fcut = 0.0_dp
      dfcutdr = 0.0_dp

      ! Cache cut_type and alias the symfgrp record locally.
      cut_type = nnp%cut_type
      ASSOCIATE (sg => nnp%rad(ind)%symfgrp(s))

      !Calculate cutoff function and partial derivative
      funccut = sg%cutoff !cutoff
      SELECT CASE (cut_type)
      CASE (nnp_cut_cos)
         tmp = r * sg%pi_over_cut
         IF (PRESENT(force)) THEN
            CALL minimax_cos_d(tmp, cos_val, dcos_dx)
            fcut = 0.5_dp*(cos_val + 1.0_dp)
            dfcutdr = 0.5_dp*dcos_dx*(pi/funccut)
         ELSE
            fcut = 0.5_dp*(minimax_cos(tmp) + 1.0_dp)
         END IF

      CASE (nnp_cut_tanh)
         tmp = 1.0_dp - r * sg%inv_cut
         tanh_val = minimax_tanh(tmp)
         fcut = tanh_val**3
         IF (PRESENT(force)) &
            dfcutdr = (-3.0_dp/funccut)*(tanh_val**2 - tanh_val**4)

      CASE DEFAULT
         CPABORT("NNP| Cutoff function unknown")
      END SELECT

      ! 1/r hoisted out of the sf loop; the three components reuse it.
      IF (PRESENT(force)) THEN
         rinv = 1.0_dp/r
         drdx(1) = rvect(1)*rinv
         drdx(2) = rvect(2)*rinv
         drdx(3) = rvect(3)*rinv
      END IF

      !loop over sym fncts of sym fnct group s
      DO sf = 1, sg%n_symf
         k = sg%symf(sf) !symf indice
         eta = nnp%rad(ind)%eta(k) !eta
         rs = nnp%rad(ind)%rs(k) !rshift

         tmp = -eta*(r - rs)**2

         ! Reuse the spline value in the chain rule (d/dx exp(x) = exp(x));
         ! more accurate than evaluating the spline derivative separately.
         IF (tmp <= -50.0_dp) THEN
            sym(sf) = 0.0_dp
            IF (PRESENT(force)) dsymdr = 0.0_dp
         ELSE IF (tmp >= 0.0_dp) THEN
            sym(sf) = 1.0_dp
            IF (PRESENT(force)) dsymdr = 0.0_dp
         ELSE
            sym(sf) = exp_reduced(tmp)
            IF (PRESENT(force)) dsymdr = sym(sf)*(-2.0_dp*eta*(r - rs))
         END IF

         ! combine radial symmetry function and cutoff function
         IF (PRESENT(force)) THEN
            force(:, sf) = fcut*dsymdr*drdx(:) + sym(sf)*dfcutdr*drdx(:)
         END IF
         sym(sf) = sym(sf)*fcut
      END DO
      END ASSOCIATE

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
   SUBROUTINE nnp_calc_ang(nnp, ind, s, rvect1, rvect2, rvect3, r1, r2, r3, sym, force, &
                            fc1_in, fc2_in, dfc1dr_in, dfc2dr_in, dr1dx_in, dr2dx_in)

      TYPE(nnp_type), INTENT(IN)                         :: nnp
      INTEGER, INTENT(IN)                                :: ind, s
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rvect1, rvect2, rvect3
      REAL(KIND=dp), INTENT(IN)                          :: r1, r2, r3
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: sym
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(OUT), &
         OPTIONAL                                        :: force
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: fc1_in, fc2_in
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: dfc1dr_in, dfc2dr_in
      REAL(KIND=dp), DIMENSION(3), INTENT(IN), OPTIONAL  :: dr1dx_in, dr2dx_in

      INTEGER                                            :: i, m, sf, cut_type
      REAL(KIND=dp) :: angular, costheta, dfcutdr1, dfcutdr2, dfcutdr3, dsymdr1, dsymdr2, dsymdr3, &
         eta, f, fcut, fcut1, fcut2, fcut3, ftot, g, g_sq_inv, lam, pref, prefzeta, rsqr1, rsqr2, &
         rsqr3, r_sum_sq, r1inv, r2inv, r3inv, symtmp, tmp, tmp1, tmp2, tmp3, tmpzeta, zeta
      REAL(KIND=dp) :: cos_val1, cos_val2, cos_val3
      REAL(KIND=dp) :: dcos_dx1, dcos_dx2, dcos_dx3
      REAL(KIND=dp) :: tanh_val1, tanh_val2, tanh_val3
      REAL(KIND=dp) :: pz2, pz4, pz8
      REAL(KIND=dp), DIMENSION(3) :: dangulardx1, dangulardx2, dangulardx3, dcosthetadx1, &
         dcosthetadx2, dcosthetadx3, dfdx1, dfdx2, dfdx3, dgdx1, dgdx2, dgdx3, dr1dx, dr2dx, dr3dx

      cut_type = nnp%cut_type
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
      ASSOCIATE (sg => nnp%ang(ind)%symfgrp(s))
      fcut = sg%cutoff !cutoff
      SELECT CASE (cut_type)
      CASE (nnp_cut_cos)
         tmp3 = r3 * sg%pi_over_cut
         IF (PRESENT(force)) THEN
            CALL minimax_cos_d(tmp3, cos_val3, dcos_dx3)
            fcut3 = 0.5_dp*(cos_val3 + 1.0_dp)
            IF (PRESENT(fc1_in)) THEN
               fcut1 = fc1_in;  fcut2 = fc2_in
            ELSE
               tmp1 = r1 * sg%pi_over_cut
               tmp2 = r2 * sg%pi_over_cut
               CALL minimax_cos_d(tmp1, cos_val1, dcos_dx1)
               CALL minimax_cos_d(tmp2, cos_val2, dcos_dx2)
               fcut1 = 0.5_dp*(cos_val1 + 1.0_dp)
               fcut2 = 0.5_dp*(cos_val2 + 1.0_dp)
            END IF
         ELSE
            fcut3 = 0.5_dp*(minimax_cos(tmp3) + 1.0_dp)
            IF (PRESENT(fc1_in)) THEN
               fcut1 = fc1_in;  fcut2 = fc2_in
            ELSE
               tmp1 = r1 * sg%pi_over_cut
               tmp2 = r2 * sg%pi_over_cut
               fcut1 = 0.5_dp*(minimax_cos(tmp1) + 1.0_dp)
               fcut2 = 0.5_dp*(minimax_cos(tmp2) + 1.0_dp)
            END IF
         END IF
         ftot = fcut1*fcut2*fcut3
         IF (PRESENT(force)) THEN
            pref = 0.5_dp*(pi/fcut)
            IF (PRESENT(fc1_in)) THEN
               dfcutdr1 = dfc1dr_in*fcut2*fcut3
               dfcutdr2 = dfc2dr_in*fcut1*fcut3
            ELSE
               dfcutdr1 = pref*dcos_dx1*fcut2*fcut3
               dfcutdr2 = pref*dcos_dx2*fcut1*fcut3
            END IF
            dfcutdr3 = pref*dcos_dx3*fcut1*fcut2
         END IF

      CASE (nnp_cut_tanh)
         tmp3 = 1.0_dp - r3 * sg%inv_cut
         tanh_val3 = minimax_tanh(tmp3)
         fcut3 = tanh_val3**3
         IF (PRESENT(fc1_in)) THEN
            fcut1 = fc1_in;  fcut2 = fc2_in
         ELSE
            tmp1 = 1.0_dp - r1 * sg%inv_cut
            tmp2 = 1.0_dp - r2 * sg%inv_cut
            tanh_val1 = minimax_tanh(tmp1)
            tanh_val2 = minimax_tanh(tmp2)
            fcut1 = tanh_val1**3
            fcut2 = tanh_val2**3
         END IF
         ftot = fcut1*fcut2*fcut3
         IF (PRESENT(force)) THEN
            pref = -3.0_dp/fcut
            IF (PRESENT(fc1_in)) THEN
               dfcutdr1 = dfc1dr_in*fcut2*fcut3
               dfcutdr2 = dfc2dr_in*fcut1*fcut3
            ELSE
               dfcutdr1 = pref*(tanh_val1**2 - tanh_val1**4)*fcut2*fcut3
               dfcutdr2 = pref*(tanh_val2**2 - tanh_val2**4)*fcut1*fcut3
            END IF
            dfcutdr3 = pref*(tanh_val3**2 - tanh_val3**4)*fcut1*fcut2
         END IF

      CASE DEFAULT
         CPABORT("NNP| Cutoff function unknown")
      END SELECT
      END ASSOCIATE

      IF (PRESENT(force)) THEN
         ! Use precomputed per-j unit vectors when available (saves 2 divisions per (j,k) pair).
         IF (PRESENT(dr1dx_in)) THEN
            dr1dx(1) = dr1dx_in(1); dr1dx(2) = dr1dx_in(2); dr1dx(3) = dr1dx_in(3)
         ELSE
            r1inv = 1.0_dp/r1
            dr1dx(1) = rvect1(1)*r1inv; dr1dx(2) = rvect1(2)*r1inv; dr1dx(3) = rvect1(3)*r1inv
         END IF
         IF (PRESENT(dr2dx_in)) THEN
            dr2dx(1) = dr2dx_in(1); dr2dx(2) = dr2dx_in(2); dr2dx(3) = dr2dx_in(3)
         ELSE
            r2inv = 1.0_dp/r2
            dr2dx(1) = rvect2(1)*r2inv; dr2dx(2) = rvect2(2)*r2inv; dr2dx(3) = rvect2(3)*r2inv
         END IF
         ! r3 unit vector is per-(j, k); always recomputed here.
         r3inv = 1.0_dp/r3
         dr3dx(1) = rvect3(1)*r3inv; dr3dx(2) = rvect3(2)*r3inv; dr3dx(3) = rvect3(3)*r3inv
      END IF

      ! Loop-invariants for the symmetry-function loop below: r_sum_sq enters
      ! the EXP argument as -eta*r_sum_sq (only eta varies per sf), and
      ! g_sq_inv enters dangular/dr.
      r_sum_sq = rsqr1 + rsqr2 + rsqr3
      g_sq_inv = 1.0_dp/(g*g)

      !loop over associated sym fncts
      ASSOCIATE (sg => nnp%ang(ind)%symfgrp(s), &
                 a  => nnp%ang(ind))
      DO sf = 1, sg%n_symf
         m = sg%symf(sf)
         lam = a%lam(m) !lambda
         zeta = a%zeta(m) !zeta
         prefzeta = a%prefzeta(m) ! 2**(1-zeta)
         eta = a%eta(m) !eta

         tmp = (1.0_dp + lam*costheta)

         IF (tmp <= 0.0_dp) THEN
            sym(sf) = 0.0_dp
            IF (PRESENT(force)) force(:, :, sf) = 0.0_dp
            CYCLE
         END IF

         ! Angular symmetry function.  zeta_int / zeta_is_int are precomputed
         ! once at init.  SELECT CASE inlines the common Behler-Parrinello
         ! exponents {1, 3, 7, 15} as plain multiplies, avoiding __powr8i4.
         IF (a%zeta_is_int(m)) THEN
            SELECT CASE (a%zeta_int(m))
            CASE (1)
               tmpzeta = 1.0_dp
            CASE (2)
               tmpzeta = tmp
            CASE (4)
               tmpzeta = tmp*tmp*tmp
            CASE (8)
               pz2 = tmp*tmp
               pz4 = pz2*pz2
               tmpzeta = pz4*pz2*tmp
            CASE (16)
               pz2 = tmp*tmp
               pz4 = pz2*pz2
               pz8 = pz4*pz4
               tmpzeta = pz8*pz4*pz2*tmp
            CASE DEFAULT
               tmpzeta = tmp**(a%zeta_int(m) - 1)
            END SELECT
         ELSE
            tmpzeta = tmp**(zeta - 1.0_dp)
         END IF
         angular = tmpzeta*tmp

         ! Spline evaluation of EXP.  The chain rule below reuses symtmp
         ! (d/dx exp(x) = exp(x)), so no second spline call is needed.
         tmp = -eta*r_sum_sq
         IF (tmp <= -50.0_dp) THEN
            symtmp = 0.0_dp
         ELSE IF (tmp >= 0.0_dp) THEN
            symtmp = 1.0_dp
         ELSE
            symtmp = exp_reduced(tmp)
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
      END ASSOCIATE

   END SUBROUTINE nnp_calc_ang

! **************************************************************************************************
!> \brief In-place heap sort of arr(1:n) in ascending order.
!>
!>        Used by nnp_calc_acsf so the downstream walks over the active set
!>        stride k in increasing order, improving prefetching and TLB
!>        locality.  O(n log n) worst case, no recursion, no allocations.
!> \param arr  integer array to sort (only arr(1:n) is touched)
!> \param n    number of valid entries
!> \date   2026-04-28
!> \author Claudio Malvino (crm98@cam.ac.uk)
! **************************************************************************************************
   PURE SUBROUTINE nnp_sort_active_atoms(arr, n)
      INTEGER, DIMENSION(:), INTENT(INOUT)               :: arr
      INTEGER, INTENT(IN)                                :: n

      INTEGER                                            :: child, heap_size, last, root, &
                                                            start, tmp_val

      IF (n < 2) RETURN

      ! Phase 1: build a max-heap by sifting each non-leaf node down.
      DO start = n/2, 1, -1
         root = start
         heap_size = n
         DO WHILE (2*root <= heap_size)
            child = 2*root
            IF (child < heap_size) THEN
               IF (arr(child) < arr(child + 1)) child = child + 1
            END IF
            IF (arr(root) >= arr(child)) EXIT
            tmp_val = arr(root); arr(root) = arr(child); arr(child) = tmp_val
            root = child
         END DO
      END DO

      ! Phase 2: swap the root into the live tail and sift the new root down.
      DO last = n, 2, -1
         tmp_val = arr(1); arr(1) = arr(last); arr(last) = tmp_val
         root = 1
         heap_size = last - 1
         DO WHILE (2*root <= heap_size)
            child = 2*root
            IF (child < heap_size) THEN
               IF (arr(child) < arr(child + 1)) child = child + 1
            END IF
            IF (arr(root) >= arr(child)) EXIT
            tmp_val = arr(root); arr(root) = arr(child); arr(child) = tmp_val
            root = child
         END DO
      END DO

   END SUBROUTINE nnp_sort_active_atoms

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

      ! The EXP/COS/TANH cutoff and Gaussian evaluations now use table-free minimax polynomials
      ! (minimax_cos / minimax_tanh / exp_reduced in splines_methods); no per-run spline tables are
      ! built here any more.

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
!> \brief Build the neighbour list for centre atom i by walking the bins
!>        around its current bin in the pre-expanded-image cell list.
!>
!>        Cost is O(k) per atom where k is the average per-bin image
!>        occupancy.  No pbc() call in the inner loop: each image carries an
!>        explicit translation, so the displacement is a direct subtraction.
!>        cell_list_cache must be populated by nnp_prepare_cell_list_cache
!>        beforehand.
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
         ! Scalar centre coordinates so the inner walk has no array temporaries.
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
               ! Skip the self-image; other periodic images of i are valid.
               IF (j == i .AND. cl%image_shift(1, img) == 0 .AND. &
                   cl%image_shift(2, img) == 0 .AND. cl%image_shift(3, img) == 0) THEN
                  img = cl%next(img)
                  CYCLE
               END IF
               drx = cx - (cl%coord_primary(1, j) + cl%image_translation(1, img))
               dry = cy - (cl%coord_primary(2, j) + cl%image_translation(2, img))
               drz = cz - (cl%coord_primary(3, j) + cl%image_translation(3, img))
               r2 = drx*drx + dry*dry + drz*drz
               ! Reject against max_cut before any SQRT.
               IF (r2 >= max_cut_sq) THEN
                  img = cl%next(img)
                  CYCLE
               END IF
               j_type = nnp%ele_ind(j)
               ! Defer SQRT until at least one group accepts this image.
               norm = -1.0_dp
               DO s = 1, nnp%rad(ind)%n_symfgrp
                  IF (j_type == nnp%rad(ind)%symfgrp(s)%ele_ind(1)) THEN
                     IF (r2 < nnp%rad(ind)%symfgrp(s)%cutoff_sq) THEN
                        IF (norm < 0.0_dp) norm = SQRT(r2)
                        n_count = neighbor%n_rad(s) + 1
                        IF (n_count > persistent_n_capacity) &
                           CPABORT("NNP| persistent neighbour buffer overflow (rad)")
                        neighbor%n_rad(s) = n_count
                        neighbor%ind_rad(n_count, s) = j
                        neighbor%dist_rad(1, n_count, s) = drx
                        neighbor%dist_rad(2, n_count, s) = dry
                        neighbor%dist_rad(3, n_count, s) = drz
                        neighbor%dist_rad(4, n_count, s) = norm
                     END IF
                  END IF
               END DO
               DO s = 1, nnp%ang(ind)%n_symfgrp
                  IF (r2 < nnp%ang(ind)%symfgrp(s)%cutoff_sq) THEN
                     IF (j_type == nnp%ang(ind)%symfgrp(s)%ele_ind(1)) THEN
                        IF (norm < 0.0_dp) norm = SQRT(r2)
                        n_count = neighbor%n_ang1(s) + 1
                        IF (n_count > persistent_n_capacity) &
                           CPABORT("NNP| persistent neighbour buffer overflow (ang1)")
                        neighbor%n_ang1(s) = n_count
                        neighbor%ind_ang1(n_count, s) = j
                        neighbor%dist_ang1(1, n_count, s) = drx
                        neighbor%dist_ang1(2, n_count, s) = dry
                        neighbor%dist_ang1(3, n_count, s) = drz
                        neighbor%dist_ang1(4, n_count, s) = norm
                     END IF
                     IF (j_type == nnp%ang(ind)%symfgrp(s)%ele_ind(2)) THEN
                        IF (norm < 0.0_dp) norm = SQRT(r2)
                        n_count = neighbor%n_ang2(s) + 1
                        IF (n_count > persistent_n_capacity) &
                           CPABORT("NNP| persistent neighbour buffer overflow (ang2)")
                        neighbor%n_ang2(s) = n_count
                        neighbor%ind_ang2(n_count, s) = j
                        neighbor%dist_ang2(1, n_count, s) = drx
                        neighbor%dist_ang2(2, n_count, s) = dry
                        neighbor%dist_ang2(3, n_count, s) = drz
                        neighbor%dist_ang2(4, n_count, s) = norm
                     END IF
                  END IF
               END DO
               img = cl%next(img)
            END DO
         END DO
         END DO
         END DO
      END ASSOCIATE

   END SUBROUTINE nnp_compute_neighbors

! **************************************************************************************************
!> \brief Lazily rebuild the module-level cell-list cache.
!>
!>        Called once per force evaluation, before the per-atom
!>        nnp_compute_neighbors loop.  The cheap path is O(N): refresh
!>        PBC-wrapped coords and check the Verlet skin.  A full O(N*(2p+1)^3)
!>        rebuild is triggered only when the maximum atomic displacement
!>        exceeds 0.5*skin or when num_atoms / cell vectors / cutoff change.
!> \param nnp ...
!> \date   2026-04-24
!> \author Claudio Malvino (crm98@cam.ac.uk)
! **************************************************************************************************
   SUBROUTINE nnp_prepare_cell_list_cache(nnp)
      TYPE(nnp_type), INTENT(INOUT), POINTER             :: nnp

      CHARACTER(len=*), PARAMETER                        :: routineN = 'nnp_prepare_cell_list_cache'
      REAL(KIND=dp), PARAMETER                           :: default_skin = 1.0_dp, &
                                                            cell_tol = 1.0e-12_dp

      INTEGER                                            :: handle, j, n_capacity, n_walk_bins
      LOGICAL                                            :: need_rebuild
      REAL(KIND=dp)                                      :: dispsq, half_skin_sq, list_cutoff, &
                                                            max_disp_sq, skin
      REAL(KIND=dp), DIMENSION(3)                        :: dr
      TYPE(cell_type), POINTER                           :: cell

      CALL timeset(routineN, handle)

      NULLIFY (cell)
      CALL nnp_env_get(nnp_env=nnp, cell=cell)

      ! Cap the skin at 10% of the cutoff so the image set stays small.
      skin = MIN(default_skin, 0.1_dp*nnp%max_cut)
      list_cutoff = nnp%max_cut + skin

      ! (Re)allocate the primary-coord buffers when num_atoms changes.
      IF (.NOT. ALLOCATED(cell_list_cache%coord_primary) .OR. &
          cell_list_cache%num_atoms /= nnp%num_atoms) THEN
         IF (ALLOCATED(cell_list_cache%coord_primary)) DEALLOCATE (cell_list_cache%coord_primary)
         IF (ALLOCATED(cell_list_cache%ref_coord_primary)) DEALLOCATE (cell_list_cache%ref_coord_primary)
         ALLOCATE (cell_list_cache%coord_primary(3, nnp%num_atoms))
         ALLOCATE (cell_list_cache%ref_coord_primary(3, nnp%num_atoms))
         cell_list_cache%initialized = .FALSE.
      END IF

      ! Refresh PBC-wrapped primary positions every step.
      DO j = 1, nnp%num_atoms
         cell_list_cache%coord_primary(:, j) = pbc(nnp%coord(:, j), cell)
      END DO

      ! Rebuild triggers.
      need_rebuild = .NOT. cell_list_cache%initialized
      IF (.NOT. need_rebuild) need_rebuild = (cell_list_cache%num_atoms /= nnp%num_atoms)
      IF (.NOT. need_rebuild) need_rebuild = (ABS(cell_list_cache%exact_cutoff - nnp%max_cut) > cell_tol)
      IF (.NOT. need_rebuild) need_rebuild = (ABS(cell_list_cache%verlet_skin - skin) > cell_tol)
      IF (.NOT. need_rebuild) need_rebuild = ANY(cell_list_cache%perd /= cell%perd)
      IF (.NOT. need_rebuild) need_rebuild = ANY(ABS(cell_list_cache%hmat - cell%hmat) > cell_tol)

      IF (.NOT. need_rebuild) THEN
         ! Verlet check: rebuild if any atom has moved more than 0.5*skin.
         ! coord_primary is already PBC-wrapped; an atom crossing a box face
         ! appears to jump by the box edge and forces a conservative rebuild,
         ! which is correct and avoids a per-atom minimum-image correction.
         half_skin_sq = (0.5_dp*skin)**2
         max_disp_sq = 0.0_dp
         DO j = 1, nnp%num_atoms
            dr(:) = cell_list_cache%coord_primary(:, j) - cell_list_cache%ref_coord_primary(:, j)
            dispsq = dr(1)**2 + dr(2)**2 + dr(3)**2
            IF (dispsq > max_disp_sq) max_disp_sq = dispsq
            IF (max_disp_sq > half_skin_sq) THEN
               need_rebuild = .TRUE.
               EXIT
            END IF
         END DO
      END IF

      IF (need_rebuild) CALL nnp_build_cell_list_cache(nnp, cell, nnp%max_cut, list_cutoff, skin)

      ! Size the persistent per-centre neighbour buffer to what one cell-list
      ! walk can actually visit (walk bins x longest chain), not to
      ! num_atoms*(2p+1)^3 images.  Chains are static between rebuilds, so
      ! the bound holds for every step the current list serves.
      n_walk_bins = MIN(cell_list_cache%nbin(1), 2*cell_list_cache%bin_span(1) + 1)* &
                    MIN(cell_list_cache%nbin(2), 2*cell_list_cache%bin_span(2) + 1)* &
                    MIN(cell_list_cache%nbin(3), 2*cell_list_cache%bin_span(3) + 1)
      n_capacity = MAX(1, MIN(cell_list_cache%n_images, &
                              n_walk_bins*cell_list_cache%max_bin_occupancy))
      CALL nnp_neighbor_ensure_persistent(nnp, n_capacity)

      CALL timestop(handle)

   END SUBROUTINE nnp_prepare_cell_list_cache

! **************************************************************************************************
!> \brief Rebuild the cell-list cache from scratch.
!>
!>        Snapshots the reference coords, expands every (atom, integer PBC
!>        shift) image, computes the padded bounding box of the image set,
!>        lays out a coarse grid with bin_width ~ list_cutoff, and links each
!>        image into its bin via the head/next arrays.
!> \param nnp ...
!> \param cell ...
!> \param exact_cutoff ...
!> \param list_cutoff ...
!> \param skin ...
!> \date   2026-04-24
!> \author Claudio Malvino (crm98@cam.ac.uk)
! **************************************************************************************************
   SUBROUTINE nnp_build_cell_list_cache(nnp, cell, exact_cutoff, list_cutoff, skin)
      TYPE(nnp_type), INTENT(IN), POINTER                :: nnp
      TYPE(cell_type), INTENT(IN), POINTER               :: cell
      REAL(KIND=dp), INTENT(IN)                          :: exact_cutoff, list_cutoff, skin

      CHARACTER(len=*), PARAMETER                        :: routineN = 'nnp_build_cell_list_cache'
      REAL(KIND=dp), PARAMETER                           :: rebuild_eps = 1.0e-10_dp

      INTEGER                                            :: cell_idx, d, handle, img, j, &
                                                            n_cells_total, n_img, n_kept, tx, ty, tz
      INTEGER, DIMENSION(3)                              :: bin
      REAL(KIND=dp), DIMENSION(3)                        :: ipos, pmax, pmin, trans, whigh, wlow

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

         ! Number of integer PBC shifts per axis that span list_cutoff.
         CALL nnp_compute_pbc_copies(cl%image_copies, cell, list_cutoff)

         ! Halo filter: centre atoms live in the primary cell, so an image
         ! can only ever be paired if its build position lies within
         ! list_cutoff of the primary bounding box (per-axis test is a
         ! conservative lower bound on the distance).  Verlet drift is
         ! covered because list_cutoff = max_cut + skin while both partners
         ! move at most skin/2 before a rebuild.  This keeps ~(1+2d/L)^3*N
         ! of the (2p+1)^3*N candidate images.
         pmin(:) = HUGE(1.0_dp)
         pmax(:) = -HUGE(1.0_dp)
         DO j = 1, nnp%num_atoms
            DO d = 1, 3
               IF (cl%coord_primary(d, j) < pmin(d)) pmin(d) = cl%coord_primary(d, j)
               IF (cl%coord_primary(d, j) > pmax(d)) pmax(d) = cl%coord_primary(d, j)
            END DO
         END DO
         wlow(:) = pmin(:) - (list_cutoff + rebuild_eps)
         whigh(:) = pmax(:) + (list_cutoff + rebuild_eps)

         ! Pass 1: count the images inside the halo window.
         n_kept = 0
         DO tz = -cl%image_copies(3), cl%image_copies(3)
         DO ty = -cl%image_copies(2), cl%image_copies(2)
         DO tx = -cl%image_copies(1), cl%image_copies(1)
            trans(:) = cell%hmat(:, 1)*REAL(tx, dp) + &
                       cell%hmat(:, 2)*REAL(ty, dp) + &
                       cell%hmat(:, 3)*REAL(tz, dp)
            DO j = 1, nnp%num_atoms
               ipos(:) = cl%coord_primary(:, j) + trans(:)
               IF (ipos(1) < wlow(1) .OR. ipos(1) > whigh(1)) CYCLE
               IF (ipos(2) < wlow(2) .OR. ipos(2) > whigh(2)) CYCLE
               IF (ipos(3) < wlow(3) .OR. ipos(3) > whigh(3)) CYCLE
               n_kept = n_kept + 1
            END DO
         END DO
         END DO
         END DO

         n_img = n_kept
         cl%n_images = n_img

         IF (ALLOCATED(cl%image_atom)) DEALLOCATE (cl%image_atom)
         IF (ALLOCATED(cl%image_shift)) DEALLOCATE (cl%image_shift)
         IF (ALLOCATED(cl%image_translation)) DEALLOCATE (cl%image_translation)
         ALLOCATE (cl%image_atom(n_img))
         ALLOCATE (cl%image_shift(3, n_img))
         ALLOCATE (cl%image_translation(3, n_img))

         ! Pass 2: fill the image table and the bounding box of kept images.
         cl%lower(:) = HUGE(1.0_dp)
         cl%upper(:) = -HUGE(1.0_dp)
         img = 0
         DO tz = -cl%image_copies(3), cl%image_copies(3)
         DO ty = -cl%image_copies(2), cl%image_copies(2)
         DO tx = -cl%image_copies(1), cl%image_copies(1)
            trans(:) = cell%hmat(:, 1)*REAL(tx, dp) + &
                       cell%hmat(:, 2)*REAL(ty, dp) + &
                       cell%hmat(:, 3)*REAL(tz, dp)
            DO j = 1, nnp%num_atoms
               ipos(:) = cl%coord_primary(:, j) + trans(:)
               IF (ipos(1) < wlow(1) .OR. ipos(1) > whigh(1)) CYCLE
               IF (ipos(2) < wlow(2) .OR. ipos(2) > whigh(2)) CYCLE
               IF (ipos(3) < wlow(3) .OR. ipos(3) > whigh(3)) CYCLE
               img = img + 1
               cl%image_atom(img) = j
               cl%image_shift(1, img) = tx
               cl%image_shift(2, img) = ty
               cl%image_shift(3, img) = tz
               cl%image_translation(:, img) = trans(:)
               DO d = 1, 3
                  IF (ipos(d) < cl%lower(d)) cl%lower(d) = ipos(d)
                  IF (ipos(d) > cl%upper(d)) cl%upper(d) = ipos(d)
               END DO
            END DO
         END DO
         END DO
         END DO

         ! Pad the bounding box by half a skin so atoms can drift before rebuild.
         cl%lower(:) = cl%lower(:) - (0.5_dp*skin + rebuild_eps)
         cl%upper(:) = cl%upper(:) + (0.5_dp*skin + rebuild_eps)

         ! Coarse grid: bin_width ~ list_cutoff so bin_span is 1 (3x3x3 walk).
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

         ! Link each image into its bin's chain (LIFO at the head).
         DO img = 1, n_img
            j = cl%image_atom(img)
            ipos(:) = cl%coord_primary(:, j) + cl%image_translation(:, img)
            CALL nnp_cell_list_bin_from_position(ipos, bin)
            cell_idx = nnp_cell_list_linear_index(bin(1), bin(2), bin(3))
            cl%next(img) = cl%head(cell_idx)
            cl%head(cell_idx) = img
         END DO

         ! Longest chain; chains are static until the next rebuild, so this
         ! exactly bounds what any 3x3x3 (or wider) walk can visit per bin.
         cl%max_bin_occupancy = 0
         DO cell_idx = 1, n_cells_total
            d = 0
            img = cl%head(cell_idx)
            DO WHILE (img > 0)
               d = d + 1
               img = cl%next(img)
            END DO
            IF (d > cl%max_bin_occupancy) cl%max_bin_occupancy = d
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

      ! Shift count for EXPLICIT image enumeration: a pair straddling a cell
      ! face pairs a primary atom with the +-1 image of its partner even when
      ! the cutoff is far below the half-box, so every periodic axis needs
      ! |T| <= 1 + cutoff/height (height = 2*proj).  The upstream per-pair
      ! pbc() search folds displacements to the minimum image and needs
      ! copies only when the cutoff exceeds the half-box; reusing that count
      ! here silently dropped every cross-face pair along any axis longer
      ! than 2*cutoff, producing step-discontinuities in the PES whenever an
      ! atom wrapped across such a face (the N>=128 instability).
      pbc_copies(1) = FLOOR(1.0_dp + cutoff/(2.0_dp*proja))
      pbc_copies(2) = FLOOR(1.0_dp + cutoff/(2.0_dp*projb))
      pbc_copies(3) = FLOOR(1.0_dp + cutoff/(2.0_dp*projc))
      ! Apply non periodic setting
      pbc_copies(:) = pbc_copies(:)*cell%perd(:)

   END SUBROUTINE nnp_compute_pbc_copies

! **************************************************************************************************
!> \brief Publish the buffer sizes every thread needs for the current system
!>        (neighbour capacity, group counts, largest symmetry-function group,
!>        y-accumulator lengths, pbc copies).
!>
!>        Runs serially once per force evaluation (from
!>        nnp_prepare_cell_list_cache).  The per-thread allocations happen in
!>        nnp_acsf_ensure_thread_buffers, called by each thread at the top of
!>        nnp_calc_acsf, because THREADPRIVATE storage can only be touched by
!>        its owning thread.
!> \param nnp ...
!> \param n_capacity  per-group row capacity; the caller derives it from the
!>                    cell-list walk bound (bins visited x longest chain)
!> \date   2026-04-28
!> \author Claudio Malvino (crm98@cam.ac.uk)
! **************************************************************************************************
   SUBROUTINE nnp_neighbor_ensure_persistent(nnp, n_capacity)
      TYPE(nnp_type), INTENT(IN), POINTER                :: nnp
      INTEGER, INTENT(IN)                                :: n_capacity

      CHARACTER(len=*), PARAMETER                        :: routineN = 'nnp_neighbor_ensure_persistent'

      INTEGER                                            :: handle, ie, n_ang_grp, &
                                 n_rad_grp, max_sym, s
      TYPE(cell_type), POINTER                           :: cell

      CALL timeset(routineN, handle)

      NULLIFY (cell)
      CALL nnp_env_get(nnp_env=nnp, cell=cell)

      CALL nnp_compute_pbc_copies(req_pbc_copies, cell, nnp%max_cut)

      n_rad_grp = 0
      n_ang_grp = 0
      DO ie = 1, nnp%n_ele
         IF (nnp%rad(ie)%n_symfgrp > n_rad_grp) n_rad_grp = nnp%rad(ie)%n_symfgrp
         IF (nnp%ang(ie)%n_symfgrp > n_ang_grp) n_ang_grp = nnp%ang(ie)%n_symfgrp
      END DO

      ! Largest symmetry-function group over all element types.
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

      req_n_capacity = MAX(1, n_capacity)
      req_n_rad_grp = MAX(1, n_rad_grp)
      req_n_ang_grp = MAX(1, n_ang_grp)
      req_max_sym = max_sym
      req_max_rad = MAX(1, MAXVAL(nnp%n_rad))
      req_max_ang = MAX(1, MAXVAL(nnp%n_ang))

      CALL timestop(handle)

   END SUBROUTINE nnp_neighbor_ensure_persistent

! **************************************************************************************************
!> \brief Grow the calling thread's persistent buffers (neighbour lists,
!>        per-pair scratch, fc precompute arrays, y accumulators) to the
!>        req_* sizes published by nnp_neighbor_ensure_persistent.
!>
!>        Fast path is a handful of integer compares; allocations happen on
!>        a thread's first centre and after system-size changes only.
!> \date   2026-06-12
!> \author Claudio Malvino (crm98@cam.ac.uk)
! **************************************************************************************************
   SUBROUTINE nnp_acsf_ensure_thread_buffers()

      IF (req_n_capacity < 1 .OR. req_max_sym < 1) &
         CPABORT("NNP| nnp_neighbor_ensure_persistent must run before nnp_calc_acsf")

      IF (scratch_max_sym < req_max_sym) THEN
         IF (ALLOCATED(scratch_sym)) &
            DEALLOCATE (scratch_sym, scratch_forcetmp, scratch_force3tmp)
         ALLOCATE (scratch_sym(req_max_sym))
         ALLOCATE (scratch_forcetmp(3, req_max_sym))
         ALLOCATE (scratch_force3tmp(3, 3, req_max_sym))
         scratch_max_sym = req_max_sym
      END IF

      IF (.NOT. ALLOCATED(ws_y_rad)) THEN
         ALLOCATE (ws_y_rad(req_max_rad), ws_y_ang(req_max_ang))
      ELSE IF (SIZE(ws_y_rad) < req_max_rad .OR. SIZE(ws_y_ang) < req_max_ang) THEN
         DEALLOCATE (ws_y_rad, ws_y_ang)
         ALLOCATE (ws_y_rad(req_max_rad), ws_y_ang(req_max_ang))
      END IF

      ! fc precompute arrays hold one entry per ang1/ang2 neighbour, which is
      ! bounded by the per-group neighbour capacity.
      IF (scratch_max_fc < req_n_capacity) THEN
         IF (ALLOCATED(scratch_fc1)) &
            DEALLOCATE (scratch_fc1, scratch_dfc1, scratch_fc2, scratch_dfc2, &
                        scratch_dr1dx, scratch_dr2dx)
         ALLOCATE (scratch_fc1(req_n_capacity), scratch_dfc1(req_n_capacity))
         ALLOCATE (scratch_fc2(req_n_capacity), scratch_dfc2(req_n_capacity))
         ALLOCATE (scratch_dr1dx(3, req_n_capacity), scratch_dr2dx(3, req_n_capacity))
         scratch_max_fc = req_n_capacity
      END IF

      IF (persistent_n_capacity < req_n_capacity .OR. &
          persistent_n_rad_grp < req_n_rad_grp .OR. &
          persistent_n_ang_grp < req_n_ang_grp) THEN
         IF (ALLOCATED(persistent_neighbor%dist_rad)) THEN
            DEALLOCATE (persistent_neighbor%dist_rad, persistent_neighbor%dist_ang1, &
                        persistent_neighbor%dist_ang2, persistent_neighbor%ind_rad, &
                        persistent_neighbor%ind_ang1, persistent_neighbor%ind_ang2, &
                        persistent_neighbor%n_rad, persistent_neighbor%n_ang1, &
                        persistent_neighbor%n_ang2)
         END IF
         ALLOCATE (persistent_neighbor%dist_rad(4, req_n_capacity, req_n_rad_grp))
         ALLOCATE (persistent_neighbor%dist_ang1(4, req_n_capacity, req_n_ang_grp))
         ALLOCATE (persistent_neighbor%dist_ang2(4, req_n_capacity, req_n_ang_grp))
         ALLOCATE (persistent_neighbor%ind_rad(req_n_capacity, req_n_rad_grp))
         ALLOCATE (persistent_neighbor%ind_ang1(req_n_capacity, req_n_ang_grp))
         ALLOCATE (persistent_neighbor%ind_ang2(req_n_capacity, req_n_ang_grp))
         ALLOCATE (persistent_neighbor%n_rad(req_n_rad_grp))
         ALLOCATE (persistent_neighbor%n_ang1(req_n_ang_grp))
         ALLOCATE (persistent_neighbor%n_ang2(req_n_ang_grp))
         persistent_n_capacity = req_n_capacity
         persistent_n_rad_grp = req_n_rad_grp
         persistent_n_ang_grp = req_n_ang_grp
      END IF

      persistent_neighbor%pbc_copies = req_pbc_copies

   END SUBROUTINE nnp_acsf_ensure_thread_buffers

! **************************************************************************************************
!> \brief Evaluate a Chebyshev series sum_k c(k) T_k(t), t in [-1,1], by Clenshaw recurrence.
!>        Value only (fast path for energy-only evaluations).
!> \param c   Chebyshev coefficients, indexed 0:n
!> \param t   evaluation point mapped to [-1, 1]
!> \return    series value at t
!> \author Claudio Malvino (crm98@cam.ac.uk)
! **************************************************************************************************
   PURE FUNCTION clenshaw(c, t) RESULT(val)
      !DIR$ ATTRIBUTES FORCEINLINE :: clenshaw
      REAL(KIND=dp), DIMENSION(0:), INTENT(IN)           :: c
      REAL(KIND=dp), INTENT(IN)                          :: t
      REAL(KIND=dp)                                      :: val

      INTEGER                                            :: j
      REAL(KIND=dp)                                      :: b0, b1, b2, t2

      t2 = 2.0_dp*t
      b1 = 0.0_dp
      b2 = 0.0_dp
      DO j = UBOUND(c, 1), 1, -1
         b0 = t2*b1 - b2 + c(j)
         b2 = b1
         b1 = b0
      END DO
      val = t*b1 - b2 + c(0)
   END FUNCTION clenshaw

! **************************************************************************************************
!> \brief As clenshaw, but also returns the derivative dval/dt of the series.  Carrying the
!>        derivative through the same recurrence guarantees the returned slope is the exact
!>        analytic derivative of the returned value, so energies and forces stay consistent.
!> \param c       Chebyshev coefficients, indexed 0:n
!> \param t       evaluation point mapped to [-1, 1]
!> \param val     series value at t
!> \param dvaldt  derivative of the series w.r.t. t at t
!> \author Claudio Malvino (crm98@cam.ac.uk)
! **************************************************************************************************
   PURE SUBROUTINE clenshaw_d(c, t, val, dvaldt)
      !DIR$ ATTRIBUTES FORCEINLINE :: clenshaw_d
      REAL(KIND=dp), DIMENSION(0:), INTENT(IN)           :: c
      REAL(KIND=dp), INTENT(IN)                          :: t
      REAL(KIND=dp), INTENT(OUT)                         :: val, dvaldt

      INTEGER                                            :: j
      REAL(KIND=dp)                                      :: b0, b1, b2, db0, db1, db2, t2

      t2 = 2.0_dp*t
      b1 = 0.0_dp; b2 = 0.0_dp
      db1 = 0.0_dp; db2 = 0.0_dp
      DO j = UBOUND(c, 1), 1, -1
         b0  = t2*b1 - b2 + c(j)
         db0 = 2.0_dp*b1 + t2*db1 - db2
         b2 = b1; b1 = b0
         db2 = db1; db1 = db0
      END DO
      val    = t*b1 - b2 + c(0)
      dvaldt = b1 + t*db1 - db2
   END SUBROUTINE clenshaw_d

! **************************************************************************************************
!> \brief Table-free approximation of cos(xx) for xx in [0, pi] via a degree-13 minimax
!>        polynomial; value only.
!> \param xx   argument in [0, pi]
!> \return     approximation of cos(xx)
!> \author Claudio Malvino (crm98@cam.ac.uk)
! **************************************************************************************************
   PURE FUNCTION minimax_cos(xx) RESULT(val)
      !DIR$ ATTRIBUTES FORCEINLINE :: minimax_cos
      REAL(KIND=dp), INTENT(IN)                          :: xx
      REAL(KIND=dp)                                      :: val

      ! map [0, pi] -> [-1, 1]
      val = clenshaw(c_cos, (2.0_dp/pi)*xx - 1.0_dp)
   END FUNCTION minimax_cos

! **************************************************************************************************
!> \brief As minimax_cos, returning value and first derivative.  A separate entry point (instead
!>        of an OPTIONAL derivative argument) keeps both variants branch-free and inlinable.
!> \param xx   argument in [0, pi]
!> \param val  approximation of cos(xx)
!> \param y1   first derivative d/dxx at xx (approximately -sin(xx))
!> \author Claudio Malvino (crm98@cam.ac.uk)
! **************************************************************************************************
   PURE SUBROUTINE minimax_cos_d(xx, val, y1)
      !DIR$ ATTRIBUTES FORCEINLINE :: minimax_cos_d
      REAL(KIND=dp), INTENT(IN)                          :: xx
      REAL(KIND=dp), INTENT(OUT)                         :: val, y1

      REAL(KIND=dp)                                      :: dvaldt

      ! map [0, pi] -> [-1, 1]
      CALL clenshaw_d(c_cos, (2.0_dp/pi)*xx - 1.0_dp, val, dvaldt)
      y1 = dvaldt*(2.0_dp/pi)
   END SUBROUTINE minimax_cos_d

! **************************************************************************************************
!> \brief Table-free approximation of tanh(xx) for xx in [0, 1] via a degree-14 minimax
!>        polynomial; value only (the cutoff derivative reuses the value algebraically).
!> \param xx   argument in [0, 1]
!> \return     approximation of tanh(xx)
!> \author Claudio Malvino (crm98@cam.ac.uk)
! **************************************************************************************************
   PURE FUNCTION minimax_tanh(xx) RESULT(val)
      !DIR$ ATTRIBUTES FORCEINLINE :: minimax_tanh
      REAL(KIND=dp), INTENT(IN)                          :: xx
      REAL(KIND=dp)                                      :: val

      ! map [0, 1] -> [-1, 1]
      val = clenshaw(c_tanh, 2.0_dp*xx - 1.0_dp)
   END FUNCTION minimax_tanh

! **************************************************************************************************
!> \brief Table-free approximation of exp(xx) via range reduction exp(x) = 2^k * exp(r),
!>        r = x - k*ln2 in [-ln2/2, ln2/2], with exp(r) from a degree-9 minimax polynomial that is
!>        accurate to full RELATIVE precision.  2^k is applied exactly with the SCALE intrinsic,
!>        so the tail decays correctly (e.g. exp(-50) is right to all digits).  d/dx exp = exp, so
!>        callers reuse the value in the chain rule; no derivative entry point is needed.
!> \param xx   argument (any real; physical ACSF inputs are <= 0)
!> \return     approximation of exp(xx)
!> \author Claudio Malvino (crm98@cam.ac.uk)
! **************************************************************************************************
   PURE FUNCTION exp_reduced(xx) RESULT(val)
      !DIR$ ATTRIBUTES FORCEINLINE :: exp_reduced
      REAL(KIND=dp), INTENT(IN)                          :: xx
      REAL(KIND=dp)                                      :: val

      INTEGER                                            :: k
      REAL(KIND=dp)                                      :: r, t

      k = NINT(xx*inv_ln2)
      r = xx - REAL(k, KIND=dp)*ln2
      ! map r in [-ln2/2, ln2/2] -> [-1, 1]
      t = r*(2.0_dp*inv_ln2)
      val = SCALE(clenshaw(c_exp, t), k)
   END FUNCTION exp_reduced

END MODULE nnp_acsf
