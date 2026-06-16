#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

# Locates the symmetrix library (https://github.com/wcwitt/symmetrix), the
# tensor-free CPU backend used by CP2K's MACE manybody potential. symmetrix
# ships two static archives (libsymmetrix + its bundled libsphericart) and is
# built against C++20; it pulls in CBLAS (provided by CP2K's BLAS) and OpenMP.
#
# Set the install location with -Dsymmetrix_ROOT=<prefix> (or the SYMMETRIX_ROOT
# environment variable); the toolchain installer exports this automatically.

include(FindPackageHandleStandardArgs)
include(cp2k_utils)

cp2k_set_default_paths(SYMMETRIX "symmetrix")

cp2k_include_dirs(SYMMETRIX "mace.hpp")
cp2k_find_libraries(SYMMETRIX "symmetrix")
# sphericart is bundled inside the symmetrix install, so look for it under the
# same root that cp2k_set_default_paths resolved for symmetrix.
set(CP2K_SYMMETRIX_SPHERICART_ROOT "${CP2K_SYMMETRIX_ROOT}")
cp2k_find_libraries(SYMMETRIX_SPHERICART "sphericart")

find_package_handle_standard_args(
  Symmetrix DEFAULT_MSG CP2K_SYMMETRIX_LINK_LIBRARIES
  CP2K_SYMMETRIX_SPHERICART_LINK_LIBRARIES CP2K_SYMMETRIX_INCLUDE_DIRS)

if(CP2K_SYMMETRIX_FOUND)
  if(NOT TARGET cp2k::symmetrix)
    add_library(cp2k::symmetrix INTERFACE IMPORTED)
  endif()
  # sphericart is parallelised with OpenMP; link the C++ OpenMP runtime so the
  # omp_* symbols resolve. CP2K always finds OpenMP (REQUIRED) before this.
  set(_symmetrix_link
      "${CP2K_SYMMETRIX_LINK_LIBRARIES};${CP2K_SYMMETRIX_SPHERICART_LINK_LIBRARIES}")
  if(TARGET OpenMP::OpenMP_CXX)
    list(APPEND _symmetrix_link OpenMP::OpenMP_CXX)
  endif()
  set_target_properties(
    cp2k::symmetrix
    PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${CP2K_SYMMETRIX_INCLUDE_DIRS}"
               INTERFACE_LINK_LIBRARIES "${_symmetrix_link}")
endif()

mark_as_advanced(CP2K_SYMMETRIX_FOUND CP2K_SYMMETRIX_INCLUDE_DIRS
                 CP2K_SYMMETRIX_LINK_LIBRARIES
                 CP2K_SYMMETRIX_SPHERICART_LINK_LIBRARIES)
