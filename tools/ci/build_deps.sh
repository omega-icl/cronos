#!/bin/bash
# Native dependencies of the CRONOS wheels (run by cibuildwheel's before-all step; usable locally too).
#   tools/ci/build_deps.sh linux|macos|none epl|gpl PREFIX
# Installs Armadillo, OpenBLAS, Eigen and Boost from the platform's packages (skipped with "none"), then builds from
# source into PREFIX the libraries the wheels bundle: SuiteSparse -- ONLY the libraries the variant links (the EPL build
# no GPL component) --, SuperLU, and SUNDIALS with KLU.  Versions are those CRONOS is validated with.
set -euo pipefail
OS=${1:?linux|macos|none} VARIANT=${2:?epl|gpl} PREFIX=${3:?install prefix}
SUITESPARSE_VERSION=7.6.1  SUPERLU_VERSION=5.2.2  SUNDIALS_VERSION=7.9.0  BOOST_VERSION=1.86.0
JOBS=$( nproc 2>/dev/null || sysctl -n hw.ncpu )
case "$VARIANT" in
  epl) SS_PROJECTS="suitesparse_config;amd;colamd;btf;klu"; SS_KLU_CHOLMOD=OFF ;;   # BSD-3-Clause, LGPL-2.1-or-later
  gpl) SS_PROJECTS="suitesparse_config;amd;camd;colamd;ccolamd;btf;klu;cholmod;umfpack;spqr"; SS_KLU_CHOLMOD=ON ;;   # + GPL-2.0-or-later
  *) echo "variant must be epl or gpl" >&2; exit 2 ;;
esac
case "$OS" in
  linux)
    dnf install -y epel-release
    dnf config-manager --set-enabled powertools
    dnf install -y armadillo-devel openblas-devel eigen3-devel
    curl -sL "https://archives.boost.io/release/${BOOST_VERSION}/source/boost_${BOOST_VERSION//./_}.tar.gz" | tar xz -C /tmp
    ( cd "/tmp/boost_${BOOST_VERSION//./_}" && ./bootstrap.sh --prefix=/usr/local > /dev/null && ./b2 install --with-system -j"$JOBS" > /dev/null ) ;;
  macos) brew install boost armadillo eigen openblas ;;
  none)  ;;
  *) echo "os must be linux, macos or none" >&2; exit 2 ;;
esac
W=$( mktemp -d ); cd "$W"
fetch(){ echo "== build_deps.sh: fetching ${1##*/}"; curl -sfL "$1" | tar xz; }
# Each step's output goes to a log, printed in full on failure (so CI shows the cause, not only "exit 1").
step(){ local name=$1; shift; echo "== build_deps.sh: $name"
  if ! "$@" > "$W/$name.log" 2>&1; then echo "build_deps.sh: $name FAILED -- its log:"; tail -n 80 "$W/$name.log"; exit 1; fi; }
# CMake 4 (current runners and manylinux images) refuses projects declaring compatibility with CMake < 3.5, as
# SuperLU 5.2.2 does; this accepts them (older CMake ignores it).
POLICY=-DCMAKE_POLICY_VERSION_MINIMUM=3.5
fetch "https://github.com/DrTimothyAldenDavis/SuiteSparse/archive/refs/tags/v${SUITESPARSE_VERSION}.tar.gz"
step suitesparse-configure cmake -S "SuiteSparse-${SUITESPARSE_VERSION}" -B ss $POLICY -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PREFIX" \
      -DSUITESPARSE_ENABLE_PROJECTS="$SS_PROJECTS" -DSUITESPARSE_USE_FORTRAN=OFF -DSUITESPARSE_USE_CUDA=OFF \
      -DSUITESPARSE_DEMOS=OFF -DBUILD_TESTING=OFF -DBUILD_STATIC_LIBS=OFF \
      -DKLU_USE_CHOLMOD=$SS_KLU_CHOLMOD   # OFF: else KLU pulls in CHOLMOD (GPL modules) for its optional ordering
step suitesparse-build cmake --build ss -j"$JOBS"
step suitesparse-install cmake --install ss
fetch "https://github.com/xiaoyeli/superlu/archive/refs/tags/v${SUPERLU_VERSION}.tar.gz"
step superlu-configure cmake -S "superlu-${SUPERLU_VERSION}" -B slu $POLICY -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PREFIX" \
      -DBUILD_SHARED_LIBS=ON -DCMAKE_POSITION_INDEPENDENT_CODE=ON -Denable_internal_blaslib=OFF -Denable_tests=OFF \
      -DXSDK_ENABLE_Fortran=OFF   # its Fortran interface is not used
step superlu-build cmake --build slu -j"$JOBS"
step superlu-install cmake --install slu
fetch "https://github.com/LLNL/sundials/releases/download/v${SUNDIALS_VERSION}/sundials-${SUNDIALS_VERSION}.tar.gz"
step sundials-configure cmake -S "sundials-${SUNDIALS_VERSION}" -B sun $POLICY -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PREFIX" \
      -DENABLE_KLU=ON -DKLU_INCLUDE_DIR="$PREFIX/include/suitesparse" -DKLU_LIBRARY_DIR="$PREFIX/lib" \
      -DBUILD_STATIC_LIBS=OFF -DEXAMPLES_ENABLE_C=OFF -DEXAMPLES_INSTALL=OFF
step sundials-build cmake --build sun -j"$JOBS"
step sundials-install cmake --install sun
echo "build_deps.sh: SuiteSparse ${SUITESPARSE_VERSION} (${SS_PROJECTS}), SuperLU ${SUPERLU_VERSION}, SUNDIALS ${SUNDIALS_VERSION} (KLU) -> ${PREFIX}"
