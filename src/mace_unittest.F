!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Unit test for the symmetrix MACE binding (mace_api). Evaluates a fixed
!>        three-atom (OHH) configuration with a reference MACE model and checks
!>        the total energy against a stored reference. The fixture is a small,
!>        purpose-built (untrained) two-element H/O MACE model, so the test pins
!>        the binding's consistency rather than physics. Mirrors nequip_unittest.F.
!> \author Claudio Malvino
! **************************************************************************************************
PROGRAM mace_unittest

   USE cp_files,                        ONLY: discover_file
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE mace_api,                        ONLY: mace_model_atomic_numbers,&
                                              mace_model_compute,&
                                              mace_model_load,&
                                              mace_model_num_elements,&
                                              mace_model_r_cut,&
                                              mace_model_release,&
                                              mace_model_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

#if defined(__MACE)
   CHARACTER(LEN=default_path_length)                 :: filename
   TYPE(mace_model_type)                              :: model
   INTEGER                                            :: nel, i, k, j, ed, iO, iH
   INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atomic_numbers
   INTEGER, PARAMETER                                 :: num_nodes = 3, num_edges = 6
   INTEGER, DIMENSION(num_nodes)                      :: node_types, num_neigh
   INTEGER, DIMENSION(num_edges)                      :: neigh_indices, neigh_types
   INTEGER, DIMENSION(2, num_nodes)                   :: nbr
   INTEGER, DIMENSION(num_nodes)                      :: zatom
   REAL(KIND=dp), DIMENSION(3*num_edges)              :: xyz
   REAL(KIND=dp), DIMENSION(num_edges)                :: r
   REAL(KIND=dp), DIMENSION(num_nodes)                :: node_energies
   REAL(KIND=dp), DIMENSION(3*num_edges)              :: node_forces
   REAL(KIND=dp), DIMENSION(3, num_nodes)             :: pos
   REAL(KIND=dp), DIMENSION(3)                        :: d
   REAL(KIND=dp)                                      :: e, ref

   ! Three-atom OHH test configuration (same as symmetrix test_mace.py).
   pos(:, 1) = [0.0_dp, -2.0_dp, 0.0_dp]   ! O
   pos(:, 2) = [1.0_dp, 0.0_dp, 0.0_dp]    ! H
   pos(:, 3) = [0.0_dp, 1.0_dp, 0.0_dp]    ! H
   zatom = [8, 1, 1]

   filename = discover_file('MACE/mace-tiny-1-8.json')
   WRITE (*, *) "Loading MACE model from: "//TRIM(filename)
   CALL mace_model_load(model, filename)
   WRITE (*, *) "Cutoff [model units]: ", mace_model_r_cut(model)

   ! Map atomic numbers to the model's 0-based element indices.
   nel = mace_model_num_elements(model)
   ALLOCATE (atomic_numbers(nel))
   CALL mace_model_atomic_numbers(model, atomic_numbers)
   iO = -1; iH = -1
   DO i = 1, nel
      IF (atomic_numbers(i) == 8) iO = i - 1
      IF (atomic_numbers(i) == 1) iH = i - 1
   END DO
   CPASSERT(iO >= 0 .AND. iH >= 0)

   ! Build the fully-connected CSR neighbour graph (all pairs within r_cut).
   ! symmetrix expects the edge vector r_i - r_j (negated w.r.t. r_j - r_i).
   nbr(:, 1) = [2, 3]
   nbr(:, 2) = [1, 3]
   nbr(:, 3) = [1, 2]
   num_neigh = [2, 2, 2]
   DO i = 1, num_nodes
      node_types(i) = MERGE(iO, iH, zatom(i) == 8)
   END DO
   ed = 0
   DO i = 1, num_nodes
      DO k = 1, 2
         ed = ed + 1
         j = nbr(k, i)
         neigh_indices(ed) = j - 1
         neigh_types(ed) = MERGE(iO, iH, zatom(j) == 8)
         d = pos(:, i) - pos(:, j)
         xyz(3*ed - 2:3*ed) = d
         r(ed) = SQRT(SUM(d*d))
      END DO
   END DO

   CALL mace_model_compute(model, num_nodes, num_edges, node_types, num_neigh, &
                           neigh_indices, neigh_types, xyz, r, &
                           node_energies, node_forces)
   e = SUM(node_energies)
   ref = 0.031956406755371_dp
   WRITE (*, *) "Total energy [model units]: ", e
   WRITE (*, *) "Reference energy          : ", ref
   WRITE (*, *) "Absolute difference       : ", ABS(e - ref)

   IF (ABS(e - ref) > 1.0e-4_dp) THEN
      CPABORT("MACE unittest failed :-(")
   END IF

   CALL mace_model_release(model)
   DEALLOCATE (atomic_numbers)

   WRITE (*, *) "MACE unittest was successful :-)"
#else
   WRITE (*, *) "MACE unittest skipped: CP2K compiled without the symmetrix MACE library."
#endif

END PROGRAM mace_unittest
