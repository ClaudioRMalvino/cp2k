!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief  Methods dealing with Neural Network potentials
!> \author Christoph Schran (christoph.schran@rub.de)
!> \date   2020-10-10
! **************************************************************************************************
MODULE nnp_force

   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                              cp_subsys_type
   USE cp_units,                        ONLY: cp_unit_from_cp2k
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE nnp_acsf,                        ONLY: nnp_calc_acsf,&
                                              nnp_prepare_cell_list_cache
   USE nnp_environment_types,           ONLY: nnp_arc_type,&
                                              nnp_env_get,&
                                              nnp_type
   USE nnp_model,                       ONLY: nnp_gradients,&
                                              nnp_predict
   USE particle_types,                  ONLY: particle_type
   USE periodic_table,                  ONLY: get_ptable_info
   USE physcon,                         ONLY: angstrom
   USE virial_types,                    ONLY: virial_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   LOGICAL, PARAMETER, PRIVATE :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'nnp_force'

   PUBLIC :: nnp_calc_energy_force

CONTAINS

! **************************************************************************************************
!> \brief Calculate the energy and force for a given configuration with the NNP
!> \param nnp ...
!> \param calc_forces ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_calc_energy_force(nnp, calc_forces)
      TYPE(nnp_type), INTENT(INOUT), POINTER             :: nnp
      LOGICAL, INTENT(IN)                                :: calc_forces

      CHARACTER(len=*), PARAMETER :: routineN = 'nnp_calc_energy_force'

      INTEGER                                            :: handle, i, i_com, idx, ie, ig, il, ind, &
                                                            istart, j, k, m, max_nsym, mecalc, &
                                                            n_active, n_base, n_remainder
      LOGICAL                                            :: calc_stress, my_expol, need_per_member
      ! Per-thread persistent state for the centre-level OMP loop: every
      ! thread owns its dsymdxyz scratch, active set, network-state clone
      ! (my_arc: node/node_grad are mutated per centre; the weight copies
      ! are read-only and cost ~100 kB/thread) and force/stress accumulators
      ! (my_force/my_cstress, folded into the shared arrays once per call).
      INTEGER, ALLOCATABLE, DIMENSION(:), SAVE                 :: active_atoms, active_flag
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), SAVE           :: denergydsym
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), SAVE    :: dsymdxyz, stress
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), SAVE    :: my_force, my_cstress
      TYPE(nnp_arc_type), ALLOCATABLE, DIMENSION(:), SAVE      :: my_arc
      !$OMP THREADPRIVATE(active_atoms, active_flag, denergydsym, dsymdxyz, stress, &
      !$OMP&              my_force, my_cstress, my_arc)
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_subsys_type), POINTER                      :: subsys
      TYPE(distribution_1d_type), POINTER                :: local_particles
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(section_vals_type), POINTER                   :: print_key, print_section
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      NULLIFY (particle_set, logger, local_particles, subsys, &
               atomic_kind_set)
      logger => cp_get_default_logger()

      CPASSERT(ASSOCIATED(nnp))
      CALL nnp_env_get(nnp_env=nnp, particle_set=particle_set, &
                       subsys=subsys, local_particles=local_particles, &
                       atomic_kind_set=atomic_kind_set)

      CALL cp_subsys_get(subsys, &
                         virial=virial)

      calc_stress = virial%pv_availability .AND. (.NOT. virial%pv_numer)
      IF (calc_stress .AND. .NOT. calc_forces) CPABORT('Stress cannot be calculated without forces')

      ! Per-committee-member forces are only consumed when the committee
      ! disagreement is used (bias) or printed; otherwise the all-reduce can
      ! carry the committee-averaged forces only, an n_committee-fold smaller
      ! message on the dominant collective.
      need_per_member = nnp%bias
      print_section => section_vals_get_subs_vals(nnp%nnp_input, "PRINT")
      print_key => section_vals_get_subs_vals(print_section, "FORCES")
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) &
         need_per_member = .TRUE.
      print_key => section_vals_get_subs_vals(print_section, "FORCES_SIGMA")
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) &
         need_per_member = .TRUE.

      ! Initialize energy and gradient:
      nnp%atomic_energy(:, :) = 0.0_dp
      IF (calc_forces) nnp%myforce(:, :, :) = 0.0_dp
      IF (calc_forces .AND. need_per_member) nnp%committee_forces(:, :, :) = 0.0_dp
      IF (calc_stress) nnp%committee_stress(:, :, :) = 0.0_dp

      ! The element-sorted index maps (sort/sort_inv/ele_ind/atoms/nuc_atoms)
      ! depend only on the particle kinds, which are fixed during MD; build
      ! them once and refresh only the coordinates afterwards.  This removes
      ! n_ele*N string comparisons and N periodic-table lookups per step.
      IF (.NOT. nnp%sort_initialized) THEN
         ig = 1
         DO i = 1, nnp%n_ele
            DO j = 1, nnp%num_atoms
               IF (nnp%ele(i) == particle_set(j)%atomic_kind%element_symbol) THEN
                  nnp%atoms(ig) = nnp%ele(i)
                  CALL get_ptable_info(nnp%atoms(ig), number=nnp%nuc_atoms(ig))
                  nnp%ele_ind(ig) = i
                  nnp%sort(ig) = j
                  nnp%sort_inv(j) = ig
                  ig = ig + 1
               END IF
            END DO
         END DO
         nnp%sort_initialized = .TRUE.
      END IF

      !fill coord array
      DO ig = 1, nnp%num_atoms
         DO m = 1, 3
            nnp%coord(m, ig) = particle_set(nnp%sort(ig))%r(m)
         END DO
      END DO

      ! Per-rank atom slice from a closed-form distribution: ranks below
      ! n_remainder get one extra atom.  Avoids the allgather present in the
      ! reference implementation.
      n_base = nnp%num_atoms/logger%para_env%num_pe
      n_remainder = MOD(nnp%num_atoms, logger%para_env%num_pe)
      IF (logger%para_env%mepos < n_remainder) THEN
         mecalc = n_base + 1
         istart = 1 + logger%para_env%mepos*(n_base + 1)
      ELSE
         mecalc = n_base
         istart = 1 + n_remainder*(n_base + 1) + (logger%para_env%mepos - n_remainder)*n_base
      END IF

      ! reset extrapolation status
      nnp%output_expol = .FALSE.

      ! Build (or refresh) the cell list once per force evaluation; the
      ! subsequent nnp_compute_neighbors calls walk only the bins around the
      ! centre atom.
      CALL nnp_prepare_cell_list_cache(nnp)

      ! dsymdxyz is dimensioned to the largest n_nodes(1) across element
      ! types so a single buffer serves every centre.
      max_nsym = 0
      DO ie = 1, nnp%n_ele
         IF (nnp%arc(ie)%n_nodes(1) > max_nsym) max_nsym = nnp%arc(ie)%n_nodes(1)
      END DO

      ! Centre-level parallelism: whole centres are distributed over the
      ! team, so nnp_calc_acsf / nnp_predict / nnp_gradients run serially
      ! within a thread on thread-owned state (no atomics, no critical
      ! sections in the hot path).  Shared writes from inside the loop are
      ! nnp%atomic_energy(i, :) only, which is disjoint per centre; force,
      ! stress, and extrapolation accumulators are per-thread and folded
      ! once at the end.
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, i_com, idx, ie, il, ind, j, k, m, n_active, my_expol)

      my_expol = .FALSE.

      ! Per-thread clone of the network state.  layer is a POINTER component,
      ! so intrinsic assignment would alias the shared layers; copy depth-2
      ! explicitly.  Weights are read-only during MD: the clone is built once
      ! and reused for the rest of the run.
      IF (.NOT. ALLOCATED(my_arc)) THEN
         ALLOCATE (my_arc(nnp%n_ele))
         DO ie = 1, nnp%n_ele
            my_arc(ie)%n_nodes = nnp%arc(ie)%n_nodes
            ALLOCATE (my_arc(ie)%layer(nnp%n_layer))
            DO il = 1, nnp%n_layer
               IF (ALLOCATED(nnp%arc(ie)%layer(il)%weights)) &
                  my_arc(ie)%layer(il)%weights = nnp%arc(ie)%layer(il)%weights
               IF (ALLOCATED(nnp%arc(ie)%layer(il)%bweights)) &
                  my_arc(ie)%layer(il)%bweights = nnp%arc(ie)%layer(il)%bweights
               IF (ALLOCATED(nnp%arc(ie)%layer(il)%node)) &
                  my_arc(ie)%layer(il)%node = nnp%arc(ie)%layer(il)%node
               IF (ALLOCATED(nnp%arc(ie)%layer(il)%node_grad)) &
                  my_arc(ie)%layer(il)%node_grad = nnp%arc(ie)%layer(il)%node_grad
            END DO
         END DO
      END IF

      IF (calc_forces) THEN
         ! Reallocate only on first call or when the system size changes.
         IF (.NOT. ALLOCATED(dsymdxyz) .OR. SIZE(dsymdxyz, 3) /= nnp%num_atoms .OR. &
             SIZE(dsymdxyz, 1) /= max_nsym) THEN
            IF (ALLOCATED(dsymdxyz)) THEN
               DEALLOCATE (dsymdxyz, active_atoms, active_flag, denergydsym)
               IF (ALLOCATED(stress)) DEALLOCATE (stress)
               IF (ALLOCATED(my_force)) DEALLOCATE (my_force)
            END IF
            ALLOCATE (dsymdxyz(max_nsym, 3, nnp%num_atoms), SOURCE=0.0_dp)
            ALLOCATE (active_atoms(nnp%num_atoms))
            ALLOCATE (active_flag(nnp%num_atoms), SOURCE=0)
            ALLOCATE (denergydsym(max_nsym))
         END IF
         IF (.NOT. ALLOCATED(my_force)) &
            ALLOCATE (my_force(3, nnp%num_atoms, nnp%n_committee))
         my_force(:, :, :) = 0.0_dp
         IF (calc_stress) THEN
            IF (.NOT. ALLOCATED(stress)) ALLOCATE (stress(3, 3, max_nsym))
            IF (.NOT. ALLOCATED(my_cstress)) ALLOCATE (my_cstress(3, 3, nnp%n_committee))
            my_cstress(:, :, :) = 0.0_dp
         END IF
      END IF

      ! calc atomic contribution to energy and force.  Dynamic schedule:
      ! centres of different elements carry different work (H vs O symmetry
      ! function counts) and the per-centre cost dwarfs the dispatch.
      !$OMP DO SCHEDULE(DYNAMIC)
      DO i = istart, istart + mecalc - 1

         ! determine index of atom type and offset
         ind = nnp%ele_ind(i)

         ! compute sym fnct values; nnp_scale_acsf fully assigns the input
         ! layer of my_arc(ind), so no reset is needed.
         IF (calc_forces) THEN
            IF (calc_stress) THEN
               stress(:, :, 1:nnp%arc(ind)%n_nodes(1)) = 0.0_dp
               CALL nnp_calc_acsf(nnp, i, dsymdxyz, stress, &
                                  active_atoms=active_atoms, active_flag=active_flag, &
                                  n_active=n_active, arc=my_arc(ind), expol=my_expol)
            ELSE
               CALL nnp_calc_acsf(nnp, i, dsymdxyz, &
                                  active_atoms=active_atoms, active_flag=active_flag, &
                                  n_active=n_active, arc=my_arc(ind), expol=my_expol)
            END IF
         ELSE
            CALL nnp_calc_acsf(nnp, i, arc=my_arc(ind), expol=my_expol)
         END IF

         DO i_com = 1, nnp%n_committee
            ! predict energy
            CALL nnp_predict(my_arc(ind), nnp, i_com)
            nnp%atomic_energy(i, i_com) = my_arc(ind)%layer(nnp%n_layer)%node(1) + &
                                          nnp%atom_energies(ind)

            ! predict forces; the inner loop is restricted to the active set
            ! built by nnp_calc_acsf, so the cost is O(n_active) per centre.
            IF (calc_forces) THEN
               ! nnp_gradients fully assigns denergydsym(1:n_nodes(1))
               CALL nnp_gradients(my_arc(ind), nnp, i_com, denergydsym)
               IF (calc_stress) THEN
                  DO j = 1, nnp%arc(ind)%n_nodes(1)
                     my_cstress(:, :, i_com) = my_cstress(:, :, i_com) - &
                                               denergydsym(j)*stress(:, :, j)
                  END DO
               END IF
               DO idx = 1, n_active
                  k = active_atoms(idx)
                  DO m = 1, 3
                     my_force(m, k, i_com) = my_force(m, k, i_com) - &
                        DOT_PRODUCT(denergydsym(1:nnp%arc(ind)%n_nodes(1)), &
                                    dsymdxyz(1:nnp%arc(ind)%n_nodes(1), m, k))
                  END DO
               END DO
            END IF
         END DO

         ! Zero only the slices of dsymdxyz that this centre touched and
         ! reset active_flag for the next centre (O(n_active) work).
         IF (calc_forces) THEN
            DO idx = 1, n_active
               k = active_atoms(idx)
               dsymdxyz(1:nnp%arc(ind)%n_nodes(1), :, k) = 0.0_dp
               active_flag(k) = 0
            END DO
         END IF

      END DO ! loop over num_atoms
      !$OMP END DO

      ! Fold the per-thread accumulators into the shared arrays.  The fold
      ! is O(3*N*n_committee) adds per thread once per force evaluation;
      ! thread order is runtime-dependent, so multi-threaded forces are
      ! bitwise-reproducible only for a fixed thread count.
      !$OMP CRITICAL (nnp_force_reduce)
      IF (calc_forces) nnp%myforce(:, :, :) = nnp%myforce(:, :, :) + my_force(:, :, :)
      IF (calc_stress) nnp%committee_stress(:, :, :) = nnp%committee_stress(:, :, :) + &
                                                       my_cstress(:, :, :)
      IF (my_expol) nnp%output_expol = .TRUE.
      !$OMP END CRITICAL (nnp_force_reduce)
      !$OMP END PARALLEL


      ! calculate energy: only the per-committee sums are consumed downstream,
      ! so reduce n_committee scalars instead of the (num_atoms, n_committee)
      ! array.  nnp%atomic_energy stays rank-local.
      nnp%committee_energy(:) = SUM(nnp%atomic_energy, 1)
      CALL logger%para_env%sum(nnp%committee_energy)
      nnp%nnp_potential_energy = SUM(nnp%committee_energy)/REAL(nnp%n_committee, dp)

      IF (calc_forces) THEN
         IF (need_per_member) THEN
            ! bring myforce to force array
            DO j = 1, nnp%num_atoms
               DO k = 1, 3
                  nnp%committee_forces(k, (nnp%sort(j)), :) = nnp%myforce(k, j, :)
               END DO
            END DO
            CALL logger%para_env%sum(nnp%committee_forces)
            nnp%nnp_forces(:, :) = SUM(nnp%committee_forces, DIM=3)/REAL(nnp%n_committee, dp)
         ELSE
            ! Average over the committee locally and reduce a single (3, N)
            ! array; committee_forces is left untouched.
            nnp%nnp_forces(:, :) = 0.0_dp
            DO i_com = 1, nnp%n_committee
               DO j = 1, nnp%num_atoms
                  nnp%nnp_forces(:, nnp%sort(j)) = nnp%nnp_forces(:, nnp%sort(j)) + &
                                                   nnp%myforce(:, j, i_com)
               END DO
            END DO
            nnp%nnp_forces(:, :) = nnp%nnp_forces(:, :)/REAL(nnp%n_committee, dp)
            CALL logger%para_env%sum(nnp%nnp_forces)
         END IF
         DO j = 1, nnp%num_atoms
            particle_set(j)%f(:) = nnp%nnp_forces(:, j)
         END DO
      END IF

      IF (calc_stress) THEN
         CALL logger%para_env%sum(nnp%committee_stress)
         virial%pv_virial = SUM(nnp%committee_stress, DIM=3)/REAL(nnp%n_committee, dp)
      END IF

      ! Bias the standard deviation of committee disagreement
      IF (nnp%bias) THEN
         CALL nnp_bias_sigma(nnp, calc_forces)
         nnp%nnp_potential_energy = nnp%nnp_potential_energy + nnp%bias_energy
         IF (calc_forces) THEN
            DO j = 1, nnp%num_atoms
               particle_set(j)%f(:) = particle_set(j)%f(:) + nnp%bias_forces(:, j)
            END DO
         END IF
         ! print properties if requested
         print_section => section_vals_get_subs_vals(nnp%nnp_input, "BIAS%PRINT")
         CALL nnp_print_bias(nnp, print_section)
      END IF

      ! print properties if requested
      print_section => section_vals_get_subs_vals(nnp%nnp_input, "PRINT")
      CALL nnp_print(nnp, print_section)

      CALL timestop(handle)

   END SUBROUTINE nnp_calc_energy_force

! **************************************************************************************************
!> \brief Calculate bias potential and force based on standard deviation of committee disagreement
!> \param nnp ...
!> \param calc_forces ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_bias_sigma(nnp, calc_forces)
      TYPE(nnp_type), INTENT(INOUT)                      :: nnp
      LOGICAL, INTENT(IN)                                :: calc_forces

      CHARACTER(len=*), PARAMETER                        :: routineN = 'nnp_bias_sigma'

      INTEGER                                            :: handle, i
      REAL(KIND=dp)                                      :: avrg, pref, sigma
      REAL(KIND=dp), DIMENSION(nnp%n_committee)          :: ecom

      CALL timeset(routineN, handle)

      ! init
      sigma = 0.0_dp
      nnp%bias_energy = 0.0_dp
      IF (calc_forces) nnp%bias_forces = 0.0_dp

      ! Subtract reference energy of each committee member, if requested.
      ! committee_energy itself stays untouched: it is the canonical
      ! globally-reduced value consumed by the PRINT%ENERGIES output.
      IF (nnp%bias_align) THEN
         ecom(:) = nnp%committee_energy(:) - nnp%bias_e_avrg(:)
      ELSE
         ecom(:) = nnp%committee_energy(:)
      END IF

      ! <E> = 1/n sum(E_i)
      ! sigma = sqrt(1/n sum((E_i - <E>)**2))
      !       = sqrt(1/n sum(dE_i**2))
      avrg = SUM(ecom)/REAL(nnp%n_committee, dp)
      DO i = 1, nnp%n_committee
         sigma = sigma + (ecom(i) - avrg)**2
      END DO
      sigma = SQRT(sigma/REAL(nnp%n_committee, dp))
      nnp%bias_sigma = sigma

      IF (sigma > nnp%bias_sigma0) THEN
         ! E_b = 0.5 * kb * (sigma - sigma_0)**2
         nnp%bias_energy = 0.5_dp*nnp%bias_kb*(sigma - nnp%bias_sigma0)**2

         IF (calc_forces) THEN
            ! nabla(E_b) = kb*(sigma - sigma_0)*nabla(sigma)
            ! nabla(sigma) = 1/sigma * 1/n sum(dE_i* nabla(dE_i))
            ! nabla(dE_i) = nabla(E_i) - nabla(<E>)
            pref = nnp%bias_kb*(1.0_dp - nnp%bias_sigma0/sigma)
            DO i = 1, nnp%n_committee
               nnp%bias_forces(:, :) = nnp%bias_forces(:, :) + &
                                       (ecom(i) - avrg)* &
                                       (nnp%committee_forces(:, :, i) - nnp%nnp_forces(:, :))
            END DO
            pref = pref/REAL(nnp%n_committee, dp)
            nnp%bias_forces(:, :) = nnp%bias_forces(:, :)*pref
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE nnp_bias_sigma

! **************************************************************************************************
!> \brief Print properties according to the requests in input file
!> \param nnp ...
!> \param print_section ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_print(nnp, print_section)
      TYPE(nnp_type), INTENT(INOUT)                      :: nnp
      TYPE(section_vals_type), INTENT(IN), POINTER       :: print_section

      INTEGER                                            :: unit_nr
      LOGICAL                                            :: explicit, file_is_new
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: print_key

      NULLIFY (logger, print_key)
      logger => cp_get_default_logger()

      print_key => section_vals_get_subs_vals(print_section, "ENERGIES")
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
         unit_nr = cp_print_key_unit_nr(logger, print_key, extension=".data", &
                                        middle_name="nnp-energies", is_new_file=file_is_new)
         IF (unit_nr > 0) CALL nnp_print_energies(nnp, unit_nr, file_is_new)
         CALL cp_print_key_finished_output(unit_nr, logger, print_key)
      END IF

      print_key => section_vals_get_subs_vals(print_section, "FORCES")
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
         IF (unit_nr > 0) CALL nnp_print_forces(nnp, print_key)
      END IF

      print_key => section_vals_get_subs_vals(print_section, "FORCES_SIGMA")
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
         unit_nr = cp_print_key_unit_nr(logger, print_key, extension=".xyz", &
                                        middle_name="nnp-forces-std")
         IF (unit_nr > 0) CALL nnp_print_force_sigma(nnp, unit_nr)
         CALL cp_print_key_finished_output(unit_nr, logger, print_key)
      END IF

      ! Output structures with extrapolation warning on any processor
      CALL logger%para_env%sum(nnp%output_expol)
      IF (nnp%output_expol) THEN
         print_key => section_vals_get_subs_vals(print_section, "EXTRAPOLATION")
         IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
            unit_nr = cp_print_key_unit_nr(logger, print_key, extension=".xyz", &
                                           middle_name="nnp-extrapolation")
            IF (unit_nr > 0) CALL nnp_print_expol(nnp, unit_nr)
            CALL cp_print_key_finished_output(unit_nr, logger, print_key)
         END IF
      END IF

      print_key => section_vals_get_subs_vals(print_section, "SUM_FORCE")

      CALL section_vals_val_get(print_section, "SUM_FORCE%ATOM_LIST", &
                                explicit=explicit)
      IF (explicit) THEN
         IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
            unit_nr = cp_print_key_unit_nr(logger, print_key, extension=".dat", &
                                           middle_name="nnp-sumforce", is_new_file=file_is_new)
            IF (unit_nr > 0) CALL nnp_print_sumforces(nnp, print_section, unit_nr, file_is_new)
            CALL cp_print_key_finished_output(unit_nr, logger, print_key)
         END IF
      END IF

   END SUBROUTINE nnp_print

! **************************************************************************************************
!> \brief Print NNP energies and standard deviation sigma
!> \param nnp ...
!> \param unit_nr ...
!> \param file_is_new ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_print_energies(nnp, unit_nr, file_is_new)
      TYPE(nnp_type), INTENT(INOUT)                      :: nnp
      INTEGER, INTENT(IN)                                :: unit_nr
      LOGICAL, INTENT(IN)                                :: file_is_new

      CHARACTER(LEN=12)                                  :: fmt_string
      INTEGER                                            :: i
      REAL(KIND=dp)                                      :: std

      IF (file_is_new) THEN
         WRITE (unit_nr, "(A1,1X,A20)", ADVANCE='no') "#", "NNP Average [a.u.],"
         WRITE (unit_nr, "(A20)", ADVANCE='no') "NNP sigma [a.u.]"
         DO i = 1, nnp%n_committee
            WRITE (unit_nr, "(A17,I3)", ADVANCE='no') "NNP", i
         END DO
         WRITE (unit_nr, "(A)") ""
      END IF

      ! nnp%atomic_energy is rank-local; committee_energy carries the
      ! globally-reduced per-member sums.
      fmt_string = "(2X,  F20.9)"
      WRITE (UNIT=fmt_string(5:6), FMT="(I2)") nnp%n_committee + 2
      std = SUM((nnp%committee_energy(:) - nnp%nnp_potential_energy)**2)
      std = std/REAL(nnp%n_committee, dp)
      std = SQRT(std)
      WRITE (unit_nr, fmt_string) nnp%nnp_potential_energy, std, nnp%committee_energy(:)

   END SUBROUTINE nnp_print_energies

! **************************************************************************************************
!> \brief Print nnp forces
!> \param nnp ...
!> \param print_key ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_print_forces(nnp, print_key)
      TYPE(nnp_type), INTENT(INOUT)                      :: nnp
      TYPE(section_vals_type), INTENT(IN), POINTER       :: print_key

      CHARACTER(len=default_path_length)                 :: fmt_string, middle_name
      INTEGER                                            :: i, j, unit_nr
      TYPE(cp_logger_type), POINTER                      :: logger

      NULLIFY (logger)
      logger => cp_get_default_logger()

      ! TODO use write_particle_coordinates from particle_methods.F
      DO i = 1, nnp%n_committee
         WRITE (fmt_string, *) i
         WRITE (middle_name, "(A,A)") "nnp-forces-", ADJUSTL(TRIM(fmt_string))
         unit_nr = cp_print_key_unit_nr(logger, print_key, extension=".xyz", &
                                        middle_name=TRIM(middle_name))
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, *) nnp%num_atoms
            WRITE (unit_nr, "(A,1X,A,A,F20.9)") "NNP forces [a.u.] of committee member", &
               ADJUSTL(TRIM(fmt_string)), "energy [a.u.]=", nnp%committee_energy(i)

            fmt_string = "(A4,1X,3F20.10)"
            DO j = 1, nnp%num_atoms
               WRITE (unit_nr, fmt_string) nnp%atoms(nnp%sort_inv(j)), nnp%committee_forces(:, j, i)
            END DO
         END IF
         CALL cp_print_key_finished_output(unit_nr, logger, print_key)
      END DO

   END SUBROUTINE nnp_print_forces

! **************************************************************************************************
!> \brief Print standard deviation sigma of NNP forces
!> \param nnp ...
!> \param unit_nr ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_print_force_sigma(nnp, unit_nr)
      TYPE(nnp_type), INTENT(INOUT)                      :: nnp
      INTEGER, INTENT(IN)                                :: unit_nr

      INTEGER                                            :: i, j
      REAL(KIND=dp), DIMENSION(3)                        :: var

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, *) nnp%num_atoms
         WRITE (unit_nr, "(A,1X,A)") "NNP sigma of forces [a.u.]"

         DO i = 1, nnp%num_atoms
            var = 0.0_dp
            DO j = 1, nnp%n_committee
               var = var + (nnp%committee_forces(:, i, j) - nnp%nnp_forces(:, i))**2
            END DO
            var = var/REAL(nnp%n_committee, dp)
            var = SQRT(var)
            WRITE (unit_nr, "(A4,1X,3F20.10)") nnp%atoms(nnp%sort_inv(i)), var
         END DO
      END IF

   END SUBROUTINE nnp_print_force_sigma

! **************************************************************************************************
!> \brief Print structures with extrapolation warning
!> \param nnp ...
!> \param unit_nr ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_print_expol(nnp, unit_nr)
      TYPE(nnp_type), INTENT(INOUT)                      :: nnp
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(len=default_path_length)                 :: fmt_string
      INTEGER                                            :: i
      REAL(KIND=dp)                                      :: mass, unit_conv
      REAL(KIND=dp), DIMENSION(3)                        :: com

      nnp%expol = nnp%expol + 1
      WRITE (unit_nr, *) nnp%num_atoms
      WRITE (unit_nr, "(A,1X,I6)") "NNP extrapolation point N =", nnp%expol

      ! move to COM of solute and wrap the box
      ! coord not needed afterwards, therefore manipulation ok
      com = 0.0_dp
      mass = 0.0_dp
      DO i = 1, nnp%num_atoms
         CALL get_ptable_info(nnp%atoms(i), amass=unit_conv)
         com(:) = com(:) + nnp%coord(:, i)*unit_conv
         mass = mass + unit_conv
      END DO
      com(:) = com(:)/mass

      DO i = 1, nnp%num_atoms
         nnp%coord(:, i) = nnp%coord(:, i) - com(:)
      END DO

      ! write out coordinates
      unit_conv = cp_unit_from_cp2k(1.0_dp, TRIM("angstrom"))
      fmt_string = "(A4,1X,3F20.10)"
      DO i = 1, nnp%num_atoms
         WRITE (unit_nr, fmt_string) &
            nnp%atoms(nnp%sort_inv(i)), &
            nnp%coord(1, nnp%sort_inv(i))*unit_conv, &
            nnp%coord(2, nnp%sort_inv(i))*unit_conv, &
            nnp%coord(3, nnp%sort_inv(i))*unit_conv
      END DO

   END SUBROUTINE nnp_print_expol

! **************************************************************************************************
!> \brief Print properties number according the requests in input file
!> \param nnp ...
!> \param print_section ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_print_bias(nnp, print_section)
      TYPE(nnp_type), INTENT(INOUT)                      :: nnp
      TYPE(section_vals_type), INTENT(IN), POINTER       :: print_section

      INTEGER                                            :: unit_nr
      LOGICAL                                            :: file_is_new
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: print_key

      NULLIFY (logger, print_key)
      logger => cp_get_default_logger()

      print_key => section_vals_get_subs_vals(print_section, "BIAS_ENERGY")
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
         unit_nr = cp_print_key_unit_nr(logger, print_key, extension=".data", &
                                        middle_name="nnp-bias-energy", is_new_file=file_is_new)
         IF (unit_nr > 0) CALL nnp_print_bias_energy(nnp, unit_nr, file_is_new)
         CALL cp_print_key_finished_output(unit_nr, logger, print_key)
      END IF

      print_key => section_vals_get_subs_vals(print_section, "BIAS_FORCES")
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
         unit_nr = cp_print_key_unit_nr(logger, print_key, extension=".xyz", &
                                        middle_name="nnp-bias-forces")
         IF (unit_nr > 0) CALL nnp_print_bias_forces(nnp, unit_nr)
         CALL cp_print_key_finished_output(unit_nr, logger, print_key)
      END IF

   END SUBROUTINE nnp_print_bias

! **************************************************************************************************
!> \brief Print NNP bias energies
!> \param nnp ...
!> \param unit_nr ...
!> \param file_is_new ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_print_bias_energy(nnp, unit_nr, file_is_new)
      TYPE(nnp_type), INTENT(INOUT)                      :: nnp
      INTEGER, INTENT(IN)                                :: unit_nr
      LOGICAL, INTENT(IN)                                :: file_is_new

      CHARACTER(len=default_path_length)                 :: fmt_string
      INTEGER                                            :: i

      IF (file_is_new) THEN
         WRITE (unit_nr, "(A1)", ADVANCE='no') "#"
         WRITE (unit_nr, "(2(2X,A19))", ADVANCE='no') "Sigma [a.u.]", "Bias energy [a.u.]"
         DO i = 1, nnp%n_committee
            IF (nnp%bias_align) THEN
               WRITE (unit_nr, "(2X,A16,I3)", ADVANCE='no') "shifted E_NNP", i
            ELSE
               WRITE (unit_nr, "(2X,A16,I3)", ADVANCE='no') "E_NNP", i
            END IF
         END DO
         WRITE (unit_nr, "(A)") ""

      END IF

      WRITE (fmt_string, "(A,I3,A)") "(2X,", nnp%n_committee + 2, "(F20.9,1X))"
      ! committee_energy is no longer overwritten by the bias routine; apply
      ! the alignment shift here when it is requested.
      IF (nnp%bias_align) THEN
         WRITE (unit_nr, fmt_string) nnp%bias_sigma, nnp%bias_energy, &
            nnp%committee_energy(:) - nnp%bias_e_avrg(:)
      ELSE
         WRITE (unit_nr, fmt_string) nnp%bias_sigma, nnp%bias_energy, nnp%committee_energy
      END IF

   END SUBROUTINE nnp_print_bias_energy

! **************************************************************************************************
!> \brief Print NNP bias forces
!> \param nnp ...
!> \param unit_nr ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_print_bias_forces(nnp, unit_nr)
      TYPE(nnp_type), INTENT(INOUT)                      :: nnp
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(len=default_path_length)                 :: fmt_string
      INTEGER                                            :: i

      ! TODO use write_particle_coordinates from particle_methods.F
      WRITE (unit_nr, *) nnp%num_atoms
      WRITE (unit_nr, "(A,F20.9)") "NNP bias forces [a.u.] for bias energy [a.u]=", nnp%bias_energy

      fmt_string = "(A4,1X,3F20.10)"
      DO i = 1, nnp%num_atoms
         WRITE (unit_nr, fmt_string) nnp%atoms(nnp%sort_inv(i)), nnp%bias_forces(:, i)
      END DO

   END SUBROUTINE nnp_print_bias_forces

! **************************************************************************************************
!> \brief Print NNP summed forces
!> \param nnp ...
!> \param print_section ...
!> \param unit_nr ...
!> \param file_is_new ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_print_sumforces(nnp, print_section, unit_nr, file_is_new)
      TYPE(nnp_type), INTENT(INOUT)                      :: nnp
      TYPE(section_vals_type), INTENT(IN), POINTER       :: print_section
      INTEGER, INTENT(IN)                                :: unit_nr
      LOGICAL, INTENT(IN)                                :: file_is_new

      CHARACTER(len=default_path_length)                 :: fmt_string
      CHARACTER(LEN=default_string_length), &
         DIMENSION(:), POINTER                           :: atomlist
      INTEGER                                            :: i, ig, j, n
      REAL(KIND=dp), DIMENSION(3)                        :: rvec

      NULLIFY (atomlist)
      IF (file_is_new) THEN
         WRITE (unit_nr, "(A)") "# Summed forces [a.u.]"
      END IF

      rvec = 0.0_dp

      ! get atoms to sum over:
      CALL section_vals_val_get(print_section, "SUM_FORCE%ATOM_LIST", &
                                c_vals=atomlist)
      IF (ASSOCIATED(atomlist)) THEN
         n = SIZE(atomlist)
         DO i = 1, nnp%num_atoms
            DO j = 1, n
               ig = nnp%sort_inv(i)
               IF (TRIM(ADJUSTL(atomlist(j))) == TRIM(ADJUSTL(nnp%atoms(ig)))) THEN
                  rvec(:) = rvec(:) + nnp%nnp_forces(:, i)
               END IF
            END DO
         END DO
      END IF

      fmt_string = "(3(F20.10,1X))"
      WRITE (unit_nr, fmt_string) rvec

   END SUBROUTINE nnp_print_sumforces

END MODULE nnp_force
