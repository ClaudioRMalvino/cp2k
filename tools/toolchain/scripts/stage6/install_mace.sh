#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

# Installs the symmetrix library (https://github.com/wcwitt/symmetrix), the
# tensor-free CPU MACE backend used by CP2K's MACE manybody potential.
# symmetrix is a CMake project with no `make install`; it produces two static
# archives (libsymmetrix.a + bundled libsphericart.a) and is header-self-
# contained under libsymmetrix/source. It needs C++20 and links CBLAS+OpenMP.
#
# NOTE: the version / tarball / checksum below are placeholders pending a
# cp2k.org download mirror; update mace_ver + mace_sha256 once the tarball is
# hosted (retrieve_package fetches from the CP2K mirror and checksum-verifies).

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

mace_ver="0.0.0"  # TODO: pin to a released symmetrix version
mace_dir="symmetrix-${mace_ver}"
mace_pkg="symmetrix-${mace_ver}.tar.gz"
mace_sha256="0000000000000000000000000000000000000000000000000000000000000000" # TODO

# shellcheck source=/dev/null
source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_mace" ] && rm "${BUILDDIR}/setup_mace"

MACE_LDFLAGS=''
MACE_LIBS=''

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_mace" in
  __INSTALL__)
    echo "==================== Installing symmetrix (MACE) ===================="
    pkg_install_dir="${INSTALLDIR}/${mace_dir}"
    install_lock_file="${pkg_install_dir}/install_successful"
    mace_root="${pkg_install_dir}"
    if verify_checksums "${install_lock_file}"; then
      echo "${mace_dir} aka symmetrix is already installed, skipping it."
    else
      retrieve_package "${mace_sha256}" "${mace_pkg}"
      echo "Installing from scratch into ${pkg_install_dir}"
      [ -d ${mace_dir} ] && rm -rf ${mace_dir}
      tar -xzf ${mace_pkg}
      cd ${mace_dir}

      # The CMake project lives in the libsymmetrix subdirectory. Kokkos is the
      # GPU/threading portability layer (default ON) and is not needed for the
      # plain CPU+MPI path; C++20 is required for std::span.
      rm -rf libsymmetrix/build
      cmake -S libsymmetrix -B libsymmetrix/build \
        -DSYMMETRIX_KOKKOS=OFF \
        -DCMAKE_CXX_STANDARD=20 \
        > libsymmetrix/cmake.log 2>&1 || tail_excerpt libsymmetrix/cmake.log

      cmake --build libsymmetrix/build -j ${NPROCS:-16} \
        > libsymmetrix/make.log 2>&1 || tail_excerpt libsymmetrix/make.log

      # no make install -- lay out a standard prefix by hand.
      [ -d ${pkg_install_dir} ] && rm -rf ${pkg_install_dir}
      mkdir -p ${pkg_install_dir}/lib ${pkg_install_dir}/include
      cp -a libsymmetrix/build/libsymmetrix.a \
        libsymmetrix/build/external/sphericart/sphericart/libsphericart.a \
        ${pkg_install_dir}/lib
      # mace.hpp pulls in its siblings by bare name, so ship the whole header
      # set (the unused *_kokkos.hpp variants are harmless).
      cp -a libsymmetrix/source/*.hpp ${pkg_install_dir}/include
      #
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage6/$(basename ${SCRIPT_NAME})"
    fi
    MACE_CFLAGS="-I'${pkg_install_dir}/include'"
    # smuggle include dirs to CXXFLAGS via DFLAGS....
    MACE_DFLAGS="-D__MACE ${MACE_CFLAGS}"
    MACE_LDFLAGS="-L'${pkg_install_dir}/lib'"
    ;;
  __DONTUSE__) ;;
  *)
    echo "==================== Linking symmetrix (MACE) to user paths ===================="
    pkg_install_dir="$with_mace"
    check_dir "${pkg_install_dir}/include"
    check_dir "${pkg_install_dir}/lib"
    MACE_CFLAGS="-I'${pkg_install_dir}/include'"
    # smuggle include dirs to CXXFLAGS via DFLAGS....
    MACE_DFLAGS="-D__MACE ${MACE_CFLAGS}"
    MACE_LDFLAGS="-L'${pkg_install_dir}/lib'"
    ;;
esac

if [ "$with_mace" != "__DONTUSE__" ]; then
  MACE_LIBS='-lsymmetrix -lsphericart -lstdc++'
  cat << EOF > "${BUILDDIR}/setup_mace"
export MACE_VER="${mace_ver}"
EOF
  if [ "$with_mace" != "__SYSTEM__" ]; then
    cat << EOF >> "${BUILDDIR}/setup_mace"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path CMAKE_PREFIX_PATH "$pkg_install_dir"
EOF
    filter_setup "${BUILDDIR}/setup_mace" "${SETUPFILE}"
  fi

  cat << EOF >> "${BUILDDIR}/setup_mace"
export MACE_DFLAGS="${MACE_DFLAGS}"
export MACE_CFLAGS="${MACE_CFLAGS}"
export MACE_LDFLAGS="${MACE_LDFLAGS}"
export MACE_LIBS="${MACE_LIBS}"
export CP_DFLAGS="\${CP_DFLAGS} ${MACE_DFLAGS}"
export CP_CFLAGS="\${CP_CFLAGS} ${MACE_CFLAGS}"
export CP_LDFLAGS="\${CP_LDFLAGS} ${MACE_LDFLAGS}"
export CP_LIBS="\${CP_LIBS} ${MACE_LIBS}"
EOF
fi

load "${BUILDDIR}/setup_mace"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "mace"
