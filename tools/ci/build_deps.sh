#!/bin/bash
# Native dependencies of the CRONOS wheels (run by cibuildwheel's before-all step; usable locally too).
#   tools/ci/build_deps.sh linux|macos|none epl|gpl PREFIX
# Installs Armadillo, OpenBLAS, Eigen and Boost from the platform's packages (skipped with "none"), then builds from
# source into PREFIX the libraries the wheels bundle: SuiteSparse -- ONLY the libraries the variant links (the EPL build
# no GPL component) --, SuperLU, and SUNDIALS with KLU.  Versions are those CRONOS is validated with.
set -euo pipefail
OS=${1:?linux|macos|none} VARIANT=${2:?epl|gpl} PREFIX=${3:?install prefix}
SUITESPARSE_VERSION=7.6.1  SUPERLU_VERSION=5.2.2  SUNDIALS_VERSION=7.9.0  BOOST_VERSION=1.86.0
OPENBLAS_WIN_VERSION=0.3.28  ARMADILLO_WIN_VERSION=14.2.2      # Windows: as in pymcpp's wheels
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
  windows)
    # As pymcpp's Windows wheels: vcpkg (x64-windows), prebuilt OpenBLAS, Armadillo's headers; ClangCL builds the module.
    # SuiteSparse comes from vcpkg -- for the EPL build suitesparse-klu only (its CHOLMOD feature is off by default).
    VCPKG=${VCPKG_INSTALLATION_ROOT:-C:/vcpkg}; VCPKG=${VCPKG//\\//}     # forward slashes (the runner sets C:\vcpkg)
    case "$VARIANT" in
      epl) VCPKG_SS="suitesparse-klu" ;;
      gpl) VCPKG_SS="suitesparse-klu suitesparse-umfpack suitesparse-spqr" ;;
    esac
    "$VCPKG/vcpkg" install --triplet x64-windows boost-interval eigen3 $VCPKG_SS
    WTAR=/c/Windows/System32/tar.exe                                   # bsdtar: reads .zip and .tar.xz
    curl -sSfL -o "$PREFIX.openblas.zip" "https://github.com/OpenMathLib/OpenBLAS/releases/download/v${OPENBLAS_WIN_VERSION}/OpenBLAS-${OPENBLAS_WIN_VERSION}-x64.zip"
    mkdir -p C:/openblas && "$WTAR" -xf "$PREFIX.openblas.zip" -C C:/openblas
    lib=$( find /c/openblas -iname 'libopenblas*.lib' | head -1 ); dll=$( find /c/openblas -iname 'libopenblas*.dll' | head -1 )
    # Canonical locations C:/openblas/libopenblas.lib and C:/openblas/bin/libopenblas.dll (where the archive may already
    # have put them: copy only what is missing -- cp refuses to copy a file onto itself)
    [ -f C:/openblas/libopenblas.lib ] || cp "$lib" C:/openblas/libopenblas.lib
    mkdir -p C:/openblas/bin; [ -f C:/openblas/bin/libopenblas.dll ] || cp "$dll" C:/openblas/bin/libopenblas.dll
    curl -sSfL -o "$PREFIX.arma.tar.xz" "https://sourceforge.net/projects/arma/files/armadillo-${ARMADILLO_WIN_VERSION}.tar.xz/download"
    mkdir -p C:/armadillo && "$WTAR" -xf "$PREFIX.arma.tar.xz" -C C:/armadillo --strip-components=1 ;;
  none)  ;;
  *) echo "os must be linux, macos, windows or none" >&2; exit 2 ;;
esac
W=$( mktemp -d ); cd "$W"
fetch(){ echo "== build_deps.sh: fetching ${1##*/}"; curl -sfL "$1" | tar xz; }
# Each step's output goes to a log, printed in full on failure (so CI shows the cause, not only "exit 1").
step(){ local name=$1; shift; echo "== build_deps.sh: $name"
  if ! "$@" > "$W/$name.log" 2>&1; then echo "build_deps.sh: $name FAILED -- its log:"; tail -n 80 "$W/$name.log"; exit 1; fi; }
# --config Release on every build and install: Visual Studio's generator (Windows) is multi-configuration and would
# build Debug -- whose runtime cannot be linked into the Release module; single-configuration generators ignore it.
# CMake 4 (current runners and manylinux images) refuses projects declaring compatibility with CMake < 3.5, as
# SuperLU 5.2.2 does; this accepts them (older CMake ignores it).
POLICY=-DCMAKE_POLICY_VERSION_MINIMUM=3.5
if [ "$OS" != windows ]; then
fetch "https://github.com/DrTimothyAldenDavis/SuiteSparse/archive/refs/tags/v${SUITESPARSE_VERSION}.tar.gz"
step suitesparse-configure cmake -S "SuiteSparse-${SUITESPARSE_VERSION}" -B ss $POLICY -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PREFIX" \
      -DSUITESPARSE_ENABLE_PROJECTS="$SS_PROJECTS" -DSUITESPARSE_USE_FORTRAN=OFF -DSUITESPARSE_USE_CUDA=OFF \
      -DSUITESPARSE_DEMOS=OFF -DBUILD_TESTING=OFF -DBUILD_STATIC_LIBS=OFF \
      -DKLU_USE_CHOLMOD=$SS_KLU_CHOLMOD   # OFF: else KLU pulls in CHOLMOD (GPL modules) for its optional ordering
step suitesparse-build cmake --build ss --config Release -j"$JOBS"
step suitesparse-install cmake --install ss --config Release
fi
fetch "https://github.com/xiaoyeli/superlu/archive/refs/tags/v${SUPERLU_VERSION}.tar.gz"
# SuperLU: its Fortran interface is not used; on Windows static (it exports no DLL symbols) and against OpenBLAS.
SLU_WIN=""
[ "$OS" = windows ] && SLU_WIN="-DBUILD_SHARED_LIBS=OFF -DTPL_BLAS_LIBRARIES=C:/openblas/libopenblas.lib -DBLAS_LIBRARIES=C:/openblas/libopenblas.lib"
step superlu-configure cmake -S "superlu-${SUPERLU_VERSION}" -B slu $POLICY -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PREFIX" \
      -DBUILD_SHARED_LIBS=ON -DCMAKE_POSITION_INDEPENDENT_CODE=ON -Denable_internal_blaslib=OFF -Denable_tests=OFF \
      -DXSDK_ENABLE_Fortran=OFF $SLU_WIN
step superlu-build cmake --build slu --config Release -j"$JOBS"
step superlu-install cmake --install slu --config Release
KLU_PREFIX=$( [ "$OS" = windows ] && echo "$VCPKG/installed/x64-windows" || echo "$PREFIX" )   # KLU: vcpkg on Windows
fetch "https://github.com/LLNL/sundials/releases/download/v${SUNDIALS_VERSION}/sundials-${SUNDIALS_VERSION}.tar.gz"
step sundials-configure cmake -S "sundials-${SUNDIALS_VERSION}" -B sun $POLICY -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PREFIX" \
      -DENABLE_KLU=ON -DKLU_INCLUDE_DIR="$KLU_PREFIX/include/suitesparse" -DKLU_LIBRARY_DIR="$KLU_PREFIX/lib" \
      -DBUILD_STATIC_LIBS=OFF -DEXAMPLES_ENABLE_C=OFF -DEXAMPLES_INSTALL=OFF
step sundials-build cmake --build sun --config Release -j"$JOBS"
step sundials-install cmake --install sun --config Release
echo "build_deps.sh ($OS): SuiteSparse $( [ "$OS" = windows ] && echo "from vcpkg (${VCPKG_SS})" || echo "${SUITESPARSE_VERSION} (${SS_PROJECTS})" ), SuperLU ${SUPERLU_VERSION}, SUNDIALS ${SUNDIALS_VERSION} (KLU) -> ${PREFIX}"
