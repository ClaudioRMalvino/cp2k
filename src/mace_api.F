!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief ISO_C_BINDING interface to the symmetrix MACE library (wcwitt/symmetrix).
!>        Mirrors the structure of torch_api.F, but the backend is the tensor-free
!>        symmetrix C++ class rather than libtorch. The model is evaluated over a
!>        local CSR neighbour graph; see mace_c_api.cpp for the C shim.
!> \author Claudio Malvino
! **************************************************************************************************
MODULE mace_api
   USE ISO_C_BINDING, ONLY: C_ASSOCIATED, &
                            C_CHAR, &
                            C_DOUBLE, &
                            C_INT, &
                            C_NULL_CHAR, &
                            C_NULL_PTR, &
                            C_PTR
   USE kinds, ONLY: dp

#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   TYPE mace_model_type
      PRIVATE
      TYPE(C_PTR)                          :: c_ptr = C_NULL_PTR
   END TYPE mace_model_type

   PUBLIC :: mace_model_type
   PUBLIC :: mace_model_load, mace_model_release
   PUBLIC :: mace_model_r_cut, mace_model_num_elements
   PUBLIC :: mace_model_atomic_numbers, mace_model_atomic_energies
   PUBLIC :: mace_model_compute

CONTAINS

! **************************************************************************************************
!> \brief Loads a MACE model from a symmetrix JSON file.
!> \param model the model handle to populate
!> \param filename path to the symmetrix JSON model file
!> \author Claudio Malvino
! **************************************************************************************************
   SUBROUTINE mace_model_load(model, filename)
      TYPE(mace_model_type), INTENT(INOUT)               :: model
      CHARACTER(LEN=*), INTENT(IN)                       :: filename

#if defined(__MACE)
      INTERFACE
         SUBROUTINE mace_c_model_load(model, filename) BIND(C, name="mace_c_model_load")
            IMPORT :: C_PTR, C_CHAR
            TYPE(C_PTR)                                  :: model
            CHARACTER(kind=C_CHAR), DIMENSION(*)         :: filename
         END SUBROUTINE mace_c_model_load
      END INTERFACE

      CPASSERT(.NOT. C_ASSOCIATED(model%c_ptr))
      CALL mace_c_model_load(model=model%c_ptr, filename=TRIM(filename)//C_NULL_CHAR)
      CPASSERT(C_ASSOCIATED(model%c_ptr))
#else
      CPABORT("CP2K compiled without the symmetrix MACE library.")
      MARK_USED(model)
      MARK_USED(filename)
#endif
   END SUBROUTINE mace_model_load

! **************************************************************************************************
!> \brief Releases a MACE model and all its resources.
!> \param model the model handle to release
!> \author Claudio Malvino
! **************************************************************************************************
   SUBROUTINE mace_model_release(model)
      TYPE(mace_model_type), INTENT(INOUT)               :: model

#if defined(__MACE)
      INTERFACE
         SUBROUTINE mace_c_model_release(model) BIND(C, name="mace_c_model_release")
            IMPORT :: C_PTR
            TYPE(C_PTR), VALUE                           :: model
         END SUBROUTINE mace_c_model_release
      END INTERFACE

      IF (C_ASSOCIATED(model%c_ptr)) THEN
         CALL mace_c_model_release(model=model%c_ptr)
         model%c_ptr = C_NULL_PTR
      END IF
#else
      CPABORT("CP2K compiled without the symmetrix MACE library.")
      MARK_USED(model)
#endif
   END SUBROUTINE mace_model_release

! **************************************************************************************************
!> \brief Returns the radial cutoff of the model (in model units).
!> \param model the model handle
!> \return the cutoff radius
!> \author Claudio Malvino
! **************************************************************************************************
   FUNCTION mace_model_r_cut(model) RESULT(r_cut)
      TYPE(mace_model_type), INTENT(IN)                  :: model
      REAL(KIND=dp)                                      :: r_cut

#if defined(__MACE)
      INTERFACE
         FUNCTION mace_c_model_r_cut(model) BIND(C, name="mace_c_model_r_cut") RESULT(r)
            IMPORT :: C_PTR, C_DOUBLE
            TYPE(C_PTR), VALUE                           :: model
            REAL(kind=C_DOUBLE)                          :: r
         END FUNCTION mace_c_model_r_cut
      END INTERFACE

      CPASSERT(C_ASSOCIATED(model%c_ptr))
      r_cut = mace_c_model_r_cut(model=model%c_ptr)
#else
      CPABORT("CP2K compiled without the symmetrix MACE library.")
      MARK_USED(model)
      r_cut = 0.0_dp
#endif
   END FUNCTION mace_model_r_cut

! **************************************************************************************************
!> \brief Returns the number of chemical elements the model was trained on.
!> \param model the model handle
!> \return the number of elements
!> \author Claudio Malvino
! **************************************************************************************************
   FUNCTION mace_model_num_elements(model) RESULT(num_elements)
      TYPE(mace_model_type), INTENT(IN)                  :: model
      INTEGER                                            :: num_elements

#if defined(__MACE)
      INTERFACE
         FUNCTION mace_c_model_num_elements(model) BIND(C, name="mace_c_model_num_elements") RESULT(n)
            IMPORT :: C_PTR, C_INT
            TYPE(C_PTR), VALUE                           :: model
            INTEGER(kind=C_INT)                          :: n
         END FUNCTION mace_c_model_num_elements
      END INTERFACE

      CPASSERT(C_ASSOCIATED(model%c_ptr))
      num_elements = INT(mace_c_model_num_elements(model=model%c_ptr))
#else
      CPABORT("CP2K compiled without the symmetrix MACE library.")
      MARK_USED(model)
      num_elements = 0
#endif
   END FUNCTION mace_model_num_elements

! **************************************************************************************************
!> \brief Copies the model's atomic numbers (one per element) into atomic_numbers.
!> \param model the model handle
!> \param atomic_numbers output, must be sized num_elements
!> \author Claudio Malvino
! **************************************************************************************************
   SUBROUTINE mace_model_atomic_numbers(model, atomic_numbers)
      TYPE(mace_model_type), INTENT(IN)                  :: model
      INTEGER, DIMENSION(:), INTENT(OUT)                 :: atomic_numbers

#if defined(__MACE)
      INTEGER(kind=C_INT), DIMENSION(SIZE(atomic_numbers)) :: tmp

      INTERFACE
         SUBROUTINE mace_c_model_atomic_numbers(model, dest) &
            BIND(C, name="mace_c_model_atomic_numbers")
            IMPORT :: C_PTR, C_INT
            TYPE(C_PTR), VALUE                           :: model
            INTEGER(kind=C_INT), DIMENSION(*)            :: dest
         END SUBROUTINE mace_c_model_atomic_numbers
      END INTERFACE

      CPASSERT(C_ASSOCIATED(model%c_ptr))
      CALL mace_c_model_atomic_numbers(model=model%c_ptr, dest=tmp)
      atomic_numbers(:) = INT(tmp)
#else
      CPABORT("CP2K compiled without the symmetrix MACE library.")
      MARK_USED(model)
      atomic_numbers(:) = 0
#endif
   END SUBROUTINE mace_model_atomic_numbers

! **************************************************************************************************
!> \brief Copies the model's per-element reference energies into atomic_energies.
!> \param model the model handle
!> \param atomic_energies output, must be sized num_elements
!> \author Claudio Malvino
! **************************************************************************************************
   SUBROUTINE mace_model_atomic_energies(model, atomic_energies)
      TYPE(mace_model_type), INTENT(IN)                  :: model
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: atomic_energies

#if defined(__MACE)
      INTERFACE
         SUBROUTINE mace_c_model_atomic_energies(model, dest) &
            BIND(C, name="mace_c_model_atomic_energies")
            IMPORT :: C_PTR, C_DOUBLE
            TYPE(C_PTR), VALUE                           :: model
            REAL(kind=C_DOUBLE), DIMENSION(*)            :: dest
         END SUBROUTINE mace_c_model_atomic_energies
      END INTERFACE

      CPASSERT(C_ASSOCIATED(model%c_ptr))
      CALL mace_c_model_atomic_energies(model=model%c_ptr, dest=atomic_energies)
#else
      CPABORT("CP2K compiled without the symmetrix MACE library.")
      MARK_USED(model)
      atomic_energies(:) = 0.0_dp
#endif
   END SUBROUTINE mace_model_atomic_energies

! **************************************************************************************************
!> \brief Evaluates the model over a local CSR neighbour graph.
!> \param model the model handle
!> \param num_nodes number of local atoms (graph nodes)
!> \param num_edges total number of directed edges = SUM(num_neigh)
!> \param node_types [num_nodes] MACE element index of each node
!> \param num_neigh [num_nodes] neighbour count of each node (CSR row sizes)
!> \param neigh_indices [num_edges] node index of each edge's neighbour
!> \param neigh_types [num_edges] MACE element index of each neighbour
!> \param xyz [3*num_edges] per-edge displacement vector
!> \param r [num_edges] per-edge distance
!> \param node_energies [num_nodes] output per-node energy
!> \param node_forces [3*num_edges] output per-edge force (scattered to atoms by caller)
!> \author Claudio Malvino
! **************************************************************************************************
   SUBROUTINE mace_model_compute(model, num_nodes, num_edges, node_types, num_neigh, &
                                 neigh_indices, neigh_types, xyz, r, &
                                 node_energies, node_forces)
      TYPE(mace_model_type), INTENT(IN)                  :: model
      INTEGER, INTENT(IN)                                :: num_nodes, num_edges
      INTEGER, DIMENSION(:), INTENT(IN)                  :: node_types, num_neigh, &
                                                            neigh_indices, neigh_types
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: xyz, r
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: node_energies, node_forces

#if defined(__MACE)
      INTEGER(kind=C_INT), DIMENSION(num_nodes)          :: node_types_c, num_neigh_c
      INTEGER(kind=C_INT), DIMENSION(num_edges)          :: neigh_indices_c, neigh_types_c

      INTERFACE
         SUBROUTINE mace_c_model_compute(model, num_nodes, num_edges, node_types, num_neigh, &
                                         neigh_indices, neigh_types, xyz, r, &
                                         node_energies, node_forces) &
            BIND(C, name="mace_c_model_compute")
            IMPORT :: C_PTR, C_INT, C_DOUBLE
            TYPE(C_PTR), VALUE                           :: model
            INTEGER(kind=C_INT), VALUE                   :: num_nodes, num_edges
            INTEGER(kind=C_INT), DIMENSION(*)            :: node_types, num_neigh, &
                                                            neigh_indices, neigh_types
            REAL(kind=C_DOUBLE), DIMENSION(*)            :: xyz, r
            REAL(kind=C_DOUBLE), DIMENSION(*)            :: node_energies, node_forces
         END SUBROUTINE mace_c_model_compute
      END INTERFACE

      CPASSERT(C_ASSOCIATED(model%c_ptr))
      ! symmetrix expects 0-based element indices; the caller (manybody_mace.F)
      ! is responsible for the kind->element-index mapping. Here we only narrow
      ! the default integers to C_INT for the ABI.
      node_types_c(:) = INT(node_types(1:num_nodes), C_INT)
      num_neigh_c(:) = INT(num_neigh(1:num_nodes), C_INT)
      neigh_indices_c(:) = INT(neigh_indices(1:num_edges), C_INT)
      neigh_types_c(:) = INT(neigh_types(1:num_edges), C_INT)

      CALL mace_c_model_compute(model=model%c_ptr, &
                                num_nodes=INT(num_nodes, C_INT), &
                                num_edges=INT(num_edges, C_INT), &
                                node_types=node_types_c, &
                                num_neigh=num_neigh_c, &
                                neigh_indices=neigh_indices_c, &
                                neigh_types=neigh_types_c, &
                                xyz=xyz, r=r, &
                                node_energies=node_energies, &
                                node_forces=node_forces)
#else
      CPABORT("CP2K compiled without the symmetrix MACE library.")
      MARK_USED(model)
      MARK_USED(num_nodes)
      MARK_USED(num_edges)
      MARK_USED(node_types)
      MARK_USED(num_neigh)
      MARK_USED(neigh_indices)
      MARK_USED(neigh_types)
      MARK_USED(xyz)
      MARK_USED(r)
      node_energies(:) = 0.0_dp
      node_forces(:) = 0.0_dp
#endif
   END SUBROUTINE mace_model_compute

END MODULE mace_api
