cmake_minimum_required(VERSION 3.18)
project(cronos LANGUAGES C CXX)

# Prefer BoostConfig.cmake with newer CMake/Boost versions.
if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()

# The CRONOS headers are validated under C++17 (the classical Make setup); C++20 is accepted.
set(CRONOS_CXX_STANDARD
    17
    CACHE STRING "C++ standard for CRONOS (17 or 20)")
set_property(CACHE CRONOS_CXX_STANDARD PROPERTY STRINGS 17 20)
set(CMAKE_CXX_STANDARD ${CRONOS_CXX_STANDARD})
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# 1. Compilation Modes
message(STATUS "--------------------------------------------------------")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified")
  set(CMAKE_BUILD_TYPE
      Release
      CACHE STRING "Choose the type of build" FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
  message(STATUS "Setting build type to 'Release'")
else()
  message(STATUS "Setting build type to 'Debug'")
endif()
message(STATUS "--------------------------------------------------------")

# Match the classical Make setup: Debug uses -O0 -g; Release uses -O2 without defining NDEBUG. CMake's default
# Release flags often include -O3 -DNDEBUG, so set the cache entries explicitly.
set(CMAKE_CXX_FLAGS_DEBUG
    "-O0 -g"
    CACHE STRING "Debug flags" FORCE)
set(CMAKE_CXX_FLAGS_RELEASE
    "-O2"
    CACHE STRING "Release flags" FORCE)

# Options
option(ENABLE_HSL "Enable HSL libraries" OFF)
option(ENABLE_EXAMPLES "Build the test drivers in test/ODESLVS and test/OCFESLV" OFF)
option(ENABLE_PYTHON "Build the cronos Python extension module" ON)
option(MC__USE_FADBAD "Use the FADBAD++ forward-AD implementation" OFF)

option(CRONOS_WITH_KLU "Sparse KLU linear solver in CVODES (SuiteSparse KLU + SUNDIALS sunlinsolklu)" ON)
option(CRONOS_WITH_UMFPACK "UMFPACK factorization in OCFESLV (SuiteSparse UMFPACK)" ON)
option(CRONOS_WITH_EIGEN "Eigen sparse QR factorization in OCFESLV" ON)
option(CRONOS_WITH_SUPERLU "SuperLU sparse factorization in OCFESLV (through arma::spsolve)" ON)
# SPQR is released under GPL-2: it is OFF by default, and a binary linking it is subject to the GPL.
option(CRONOS_WITH_SPQR "SPQR factorization in OCFESLV (SuiteSparse SPQR -- GPL-2 licensed)" OFF)

set(MC_INTERVAL_LIBRARY
    "BOOST"
    CACHE STRING
          "Interval library backend: BOOST, PROFIL, FILIB, or NONVERIFIED")
set_property(CACHE MC_INTERVAL_LIBRARY PROPERTY STRINGS BOOST PROFIL FILIB
                                                NONVERIFIED)

# MC++ (required)
#
# MC++ v5 installs its headers but exports no CMake package, so CRONOS reads them from the MC++ source tree:
# MCPP_ROOT is the directory holding src/mc (cache variable or environment variable).
set(MCPP_ROOT
    "$ENV{MCPP_ROOT}"
    CACHE PATH "MC++ source tree (the directory holding src/mc)")
if(NOT MCPP_ROOT OR NOT EXISTS "${MCPP_ROOT}/src/mc/mcfunc.hpp")
  message(
    FATAL_ERROR
      "MC++ not found: set MCPP_ROOT (cache or environment variable) to the MC++ source tree, the directory "
      "holding src/mc/mcfunc.hpp (MCPP_ROOT='${MCPP_ROOT}')")
endif()
set(MCPP_INCLUDE_DIRS "${MCPP_ROOT}/src/mc")
if(MC__USE_FADBAD)
  if(NOT EXISTS "${MCPP_ROOT}/src/3rdparty/fadbad++")
    message(FATAL_ERROR "MC__USE_FADBAD=ON but ${MCPP_ROOT}/src/3rdparty/fadbad++ does not exist")
  endif()
  list(APPEND MCPP_INCLUDE_DIRS "${MCPP_ROOT}/src/3rdparty/fadbad++")
endif()
message(STATUS "MC++ headers: ${MCPP_INCLUDE_DIRS}")

# Dependencies

# Python 3 Configuration CUSTOM_PYTHON_PATH may be either: * an executable, e.g. /path/to/.venv/bin/python * a
# virtual-environment/root directory, e.g. /path/to/.venv Leave empty to let CMake find Python.
if(ENABLE_PYTHON)
  set(CUSTOM_PYTHON_PATH
      ""
      CACHE STRING "Custom Python root directory or Python executable")

  if(CUSTOM_PYTHON_PATH)
    if(EXISTS "${CUSTOM_PYTHON_PATH}" AND NOT IS_DIRECTORY
                                          "${CUSTOM_PYTHON_PATH}")
      set(Python3_EXECUTABLE
          "${CUSTOM_PYTHON_PATH}"
          CACHE FILEPATH "Python executable" FORCE)
      message(
        STATUS
          "Custom Python path recognized as executable: ${Python3_EXECUTABLE}")
    elseif(IS_DIRECTORY "${CUSTOM_PYTHON_PATH}")
      set(Python3_ROOT_DIR
          "${CUSTOM_PYTHON_PATH}"
          CACHE PATH "Python root directory" FORCE)
      find_program(
        Python3_EXECUTABLE_HINT
        NAMES python3 python
        PATHS "${CUSTOM_PYTHON_PATH}/bin" "${CUSTOM_PYTHON_PATH}"
        NO_DEFAULT_PATH NO_CMAKE_ENVIRONMENT_PATH NO_SYSTEM_ENVIRONMENT_PATH)
      if(Python3_EXECUTABLE_HINT)
        set(Python3_EXECUTABLE
            "${Python3_EXECUTABLE_HINT}"
            CACHE FILEPATH "Python executable" FORCE)
        message(
          STATUS "Found Python executable in custom root: ${Python3_EXECUTABLE}")
      else()
        message(
          FATAL_ERROR
            "CUSTOM_PYTHON_PATH='${CUSTOM_PYTHON_PATH}' is a directory, but no python/python3 executable was found in it"
        )
      endif()
    else()
      message(
        FATAL_ERROR "CUSTOM_PYTHON_PATH='${CUSTOM_PYTHON_PATH}' does not exist")
    endif()
  endif()

  # Find Python with the selected executable locked in. Development.Module is the appropriate component for building
  # extension modules; it avoids unnecessarily requesting the embedding library.
  set(Python3_FIND_VIRTUALENV FIRST)
  find_package(
    Python3
    COMPONENTS Interpreter Development.Module
    REQUIRED)

  execute_process(
    COMMAND "${Python3_EXECUTABLE}" -c
            "import sysconfig; print(sysconfig.get_path('platlib') or '')"
    OUTPUT_VARIABLE Python3_SITEARCH
    OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)

  message(STATUS "--------------------------------------------------------")
  message(STATUS "Python3 Interpreter: ${Python3_EXECUTABLE}")
  message(STATUS "Python3 Module libs: ${Python3_LIBRARIES}")
  message(STATUS "Python3 Includes:    ${Python3_INCLUDE_DIRS}")
  message(STATUS "Python3 site-packages: ${Python3_SITEARCH}")
  message(STATUS "--------------------------------------------------------")

  # pybind11 (>= 3.0.3, NOT 3.1.0; 3.0.4 recommended).  3.0.3 brought the virtual-base detection and the
  # string-ownership fix (gh-5976) CRONOS relies on; 3.1.0 has a keep_alive<0,N> regression on rejected overloads.
  # cronos uses pymcpp's types (FFVar, FFGraph, ...) across the module boundary, which pybind11 shares only between
  # modules built against the SAME internals: build pymcpp with the same pybind11 release.
  # Order: the vendored extern/pybind11, else an installed pybind11 package, else (CRONOS_FETCH_PYBIND11) v3.0.4
  # fetched from GitHub.
  option(CRONOS_FETCH_PYBIND11 "Fetch pybind11 v3.0.4 when no acceptable pybind11 is vendored or installed" ON)
  set(CRONOS_PYBIND11_TAG "v3.0.4" CACHE STRING "pybind11 release fetched by CRONOS_FETCH_PYBIND11")

  # cronos_pybind11_version(<include dir> <out var>): the release, read from pybind11/detail/common.h
  function(cronos_pybind11_version incdir outvar)
    set(_h "${incdir}/pybind11/detail/common.h")
    if(NOT EXISTS "${_h}")
      set(${outvar} "" PARENT_SCOPE)
      return()
    endif()
    file(STRINGS "${_h}" _l REGEX "#define PYBIND11_VERSION_(MAJOR|MINOR|PATCH) ")
    set(_v "")
    foreach(_k MAJOR MINOR PATCH)
      string(REGEX MATCH "PYBIND11_VERSION_${_k} ([0-9]+)" _m "${_l}")
      list(APPEND _v "${CMAKE_MATCH_1}")
    endforeach()
    string(REPLACE ";" "." _v "${_v}")
    set(${outvar} "${_v}" PARENT_SCOPE)
  endfunction()
  # cronos_pybind11_ok(<version> <out var>): >= 3.0.3 and not 3.1.0
  function(cronos_pybind11_ok ver outvar)
    if(ver AND NOT ver VERSION_LESS "3.0.3" AND NOT ver VERSION_EQUAL "3.1.0")
      set(${outvar} ON PARENT_SCOPE)
    else()
      set(${outvar} OFF PARENT_SCOPE)
    endif()
  endfunction()

  set(PYBIND11_FINDPYTHON ON)
  set(CRONOS_VENDORED_PYBIND11 OFF)
  set(_pyb_src "")
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/extern/pybind11/CMakeLists.txt")
    cronos_pybind11_version("${CMAKE_CURRENT_SOURCE_DIR}/extern/pybind11/include" _pyb_ver)
    cronos_pybind11_ok("${_pyb_ver}" _pyb_ok)
    if(NOT _pyb_ok)
      message(FATAL_ERROR
        "extern/pybind11 is pybind11 ${_pyb_ver}: CRONOS needs >= 3.0.3 and not 3.1.0 (keep_alive regression on "
        "rejected overloads).  Check out v3.0.4 there (git -C extern/pybind11 checkout v3.0.4), or remove "
        "extern/pybind11 to use an installed or fetched pybind11.")
    endif()
    add_subdirectory(extern/pybind11)
    set(CRONOS_VENDORED_PYBIND11 ON)
    set(_pyb_src "vendored (extern/pybind11)")
  else()
    find_package(pybind11 CONFIG QUIET)
    set(_pyb_ok OFF)
    if(pybind11_FOUND)
      set(_pyb_ver "${pybind11_VERSION}")
      cronos_pybind11_ok("${_pyb_ver}" _pyb_ok)
      if(NOT _pyb_ok)
        message(STATUS "Installed pybind11 ${_pyb_ver} (${pybind11_DIR}) is not acceptable (needs >= 3.0.3, not 3.1.0)")
      endif()
    endif()
    if(_pyb_ok)
      set(_pyb_src "installed (${pybind11_DIR})")
    elseif(CRONOS_FETCH_PYBIND11)
      include(FetchContent)
      FetchContent_Declare(
        pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11.git
        GIT_TAG ${CRONOS_PYBIND11_TAG}
        GIT_SHALLOW TRUE)
      FetchContent_MakeAvailable(pybind11)
      cronos_pybind11_version("${pybind11_SOURCE_DIR}/include" _pyb_ver)
      set(_pyb_src "fetched (${CRONOS_PYBIND11_TAG})")
    else()
      message(FATAL_ERROR
        "No acceptable pybind11 (>= 3.0.3, not 3.1.0): vendor v3.0.4 in extern/pybind11, install it, or pass "
        "-DCRONOS_FETCH_PYBIND11=ON")
    endif()
  endif()
  message(STATUS "pybind11 ${_pyb_ver}: ${_pyb_src} -- must be the release pymcpp was built with")

  # A changed pybind11 release is not seen by CMake's dependency tracking when it sits under the same path (a
  # retargeted extern/pybind11, an upgraded package): some units would be rebuilt against the new release and others
  # kept from the old one, and a module mixing pybind11 internals misbehaves (garbage enum values, crashes).  Refuse
  # to continue in a build directory configured with another release.
  if(DEFINED CRONOS_PYBIND11_CONFIGURED AND NOT CRONOS_PYBIND11_CONFIGURED STREQUAL "${_pyb_ver}")
    message(FATAL_ERROR
      "pybind11 changed from ${CRONOS_PYBIND11_CONFIGURED} to ${_pyb_ver} since this build directory was configured. "
      "CMake does not track that change, so start from a CLEAN build directory (and rebuild pymcpp with the same "
      "release).")
  endif()
  set(CRONOS_PYBIND11_CONFIGURED "${_pyb_ver}" CACHE INTERNAL "pybind11 release this build directory was configured with")
endif()

# Boost Always find Boost because it might be used for other things or by default
find_package(Boost CONFIG REQUIRED)

# Armadillo
find_package(Armadillo REQUIRED)

# LAPACK / BLAS
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# Threads (MC__USE_THREAD)
find_package(Threads REQUIRED)

# HSL
if(ENABLE_HSL)
  message(STATUS "Enabling HSL support")
  find_library(LIB_MC13 NAMES mc13)
  find_library(LIB_MC21 NAMES mc21)
  find_library(LIB_MC33 NAMES mc33)

  set(HSL_LIBRARIES "")
  if(LIB_MC13)
    message(STATUS "  Found MC13: ${LIB_MC13}")
    list(APPEND HSL_LIBRARIES ${LIB_MC13})
  else()
    message(
      STATUS
        "  MC13 not found via find_library, assuming 'mc13' is in linker path")
    list(APPEND HSL_LIBRARIES mc13)
  endif()

  if(LIB_MC21)
    message(STATUS "  Found MC21: ${LIB_MC21}")
    list(APPEND HSL_LIBRARIES ${LIB_MC21})
  else()
    message(
      STATUS
        "  MC21 not found via find_library, assuming 'mc21' is in linker path")
    list(APPEND HSL_LIBRARIES mc21)
  endif()

  if(LIB_MC33)
    message(STATUS "  Found MC33: ${LIB_MC33}")
    list(APPEND HSL_LIBRARIES ${LIB_MC33})
  else()
    message(
      STATUS
        "  MC33 not found via find_library, assuming 'mc33' is in linker path")
    list(APPEND HSL_LIBRARIES mc33)
  endif()

  list(APPEND HSL_LIBRARIES gfortran)
endif()

# Interval Library Configuration (as MC++)
message(STATUS "--------------------------------------------------------")
message(STATUS "Active Interval Library Backend: ${MC_INTERVAL_LIBRARY}")
message(STATUS "--------------------------------------------------------")

set(INTERVAL_LIBRARIES "")
set(INTERVAL_INCLUDE_DIRS "")
set(INTERVAL_DEFINITIONS "")
set(INTERVAL_COMPILE_OPTIONS "")

if(MC_INTERVAL_LIBRARY STREQUAL "BOOST")
  list(APPEND INTERVAL_DEFINITIONS MC__USE_BOOST)
elseif(MC_INTERVAL_LIBRARY STREQUAL "PROFIL")
  if(NOT DEFINED PROFIL_HOME)
    if(DEFINED ENV{PROFIL_HOME})
      set(PROFIL_HOME $ENV{PROFIL_HOME})
    else()
      message(
        FATAL_ERROR
          "PROFIL_HOME not set. Please set PROFIL_HOME environment variable or cmake variable."
      )
    endif()
  endif()
  list(APPEND INTERVAL_INCLUDE_DIRS ${PROFIL_HOME}/include)
  list(
    APPEND
    INTERVAL_LIBRARIES
    -L${PROFIL_HOME}/lib
    ProfilPackages
    Profil
    Bias
    lr)
  list(APPEND INTERVAL_DEFINITIONS MC__USE_PROFIL)
elseif(MC_INTERVAL_LIBRARY STREQUAL "FILIB")
  if(NOT DEFINED FILIB_HOME)
    if(DEFINED ENV{FILIB_HOME})
      set(FILIB_HOME $ENV{FILIB_HOME})
    else()
      message(
        FATAL_ERROR
          "FILIB_HOME not set. Please set FILIB_HOME environment variable or cmake variable."
      )
    endif()
  endif()
  list(APPEND INTERVAL_INCLUDE_DIRS ${FILIB_HOME}/include
       ${FILIB_HOME}/include/interval)
  list(APPEND INTERVAL_LIBRARIES -L${FILIB_HOME}/lib prim)
  list(APPEND INTERVAL_DEFINITIONS MC__USE_FILIB)
  list(APPEND INTERVAL_COMPILE_OPTIONS -frounding-math)
elseif(MC_INTERVAL_LIBRARY STREQUAL "NONVERIFIED")
  # No interval library
else()
  message(
    FATAL_ERROR
      "Invalid MC_INTERVAL_LIBRARY: ${MC_INTERVAL_LIBRARY}. Must be one of BOOST, PROFIL, FILIB, NONVERIFIED.")
endif()

# SUNDIALS (required, version >= 7: CVODES)
#
# Found WITHOUT any hint when installed in a standard place (/usr, /usr/local, CMAKE_INSTALL_PREFIX), on
# CMAKE_PREFIX_PATH, or under /opt/sundials*; otherwise point SUNDIALS_ROOT (cache or environment variable) at its
# installation prefix, or SUNDIALS_DIR at <prefix>/lib/cmake/sundials.  SUNDIALS >= 6 installs SUNDIALSConfig.cmake;
# a SUNDIALS without it is found library by library, its version read from sundials_config.h.
set(CRONOS_SUNDIALS_MIN 7)
set(SUNDIALS_ROOT
    "$ENV{SUNDIALS_ROOT}"
    CACHE PATH "Optional SUNDIALS installation prefix")
file(GLOB _sundials_opt LIST_DIRECTORIES true "/opt/sundials*" "/usr/local/sundials*")
set(_cronos_sundials_components
    cvodes
    nvecserial
    sunmatrixdense
    sunlinsoldense
    sunmatrixsparse
    sunnonlinsolnewton
    sunnonlinsolfixedpoint)
if(CRONOS_WITH_KLU)
  list(APPEND _cronos_sundials_components sunlinsolklu)
endif()

set(SUNDIALS_TARGETS "")
find_package(SUNDIALS ${CRONOS_SUNDIALS_MIN} CONFIG QUIET HINTS "${SUNDIALS_ROOT}" PATHS ${_sundials_opt})
if(NOT SUNDIALS_FOUND)
  # an older SUNDIALS on the path: say so, rather than "not found"
  find_package(SUNDIALS CONFIG QUIET HINTS "${SUNDIALS_ROOT}" PATHS ${_sundials_opt})
  if(SUNDIALS_FOUND)
    message(FATAL_ERROR
      "SUNDIALS ${SUNDIALS_VERSION} found (${SUNDIALS_DIR}), but CRONOS requires SUNDIALS >= ${CRONOS_SUNDIALS_MIN}: "
      "install SUNDIALS ${CRONOS_SUNDIALS_MIN}.x (with -DENABLE_KLU=ON for CRONOS_WITH_KLU) and point SUNDIALS_ROOT "
      "at it")
  endif()
endif()
if(SUNDIALS_FOUND)
  message(STATUS "SUNDIALS ${SUNDIALS_VERSION}: ${SUNDIALS_DIR}")
  foreach(_c IN LISTS _cronos_sundials_components)
    if(TARGET SUNDIALS::${_c})
      list(APPEND SUNDIALS_TARGETS SUNDIALS::${_c})
    elseif(TARGET SUNDIALS::${_c}_shared)
      list(APPEND SUNDIALS_TARGETS SUNDIALS::${_c}_shared)
    elseif(TARGET SUNDIALS::${_c}_static)
      list(APPEND SUNDIALS_TARGETS SUNDIALS::${_c}_static)
    else()
      message(FATAL_ERROR
        "SUNDIALS was found but provides no target for component '${_c}' -- with CRONOS_WITH_KLU=ON, SUNDIALS "
        "must be built with KLU (-DENABLE_KLU=ON)")
    endif()
  endforeach()
  foreach(_core SUNDIALS::core SUNDIALS::core_shared SUNDIALS::core_static)
    if(TARGET ${_core})
      list(APPEND SUNDIALS_TARGETS ${_core})
      break()
    endif()
  endforeach()
else()
  find_path(SUNDIALS_INCLUDE_DIR cvodes/cvodes.h HINTS "${SUNDIALS_ROOT}/include" PATHS ${_sundials_opt}
            PATH_SUFFIXES include)
  if(NOT SUNDIALS_INCLUDE_DIR)
    message(FATAL_ERROR
      "SUNDIALS (>= ${CRONOS_SUNDIALS_MIN}) not found: set SUNDIALS_ROOT to its installation prefix (cache or "
      "environment variable), or add the prefix to CMAKE_PREFIX_PATH")
  endif()
  set(_sv "")
  if(EXISTS "${SUNDIALS_INCLUDE_DIR}/sundials/sundials_config.h")
    file(STRINGS "${SUNDIALS_INCLUDE_DIR}/sundials/sundials_config.h" _svl REGEX "#define SUNDIALS_VERSION_MAJOR ")
    string(REGEX MATCH "SUNDIALS_VERSION_MAJOR ([0-9]+)" _m "${_svl}")
    set(_sv "${CMAKE_MATCH_1}")
  endif()
  if(NOT _sv OR _sv LESS CRONOS_SUNDIALS_MIN)
    message(FATAL_ERROR
      "SUNDIALS at ${SUNDIALS_INCLUDE_DIR} is version '${_sv}': CRONOS requires >= ${CRONOS_SUNDIALS_MIN}")
  endif()
  get_filename_component(_sroot "${SUNDIALS_INCLUDE_DIR}" DIRECTORY)
  foreach(_c IN LISTS _cronos_sundials_components ITEMS core)
    find_library(SUNDIALS_LIB_${_c} NAMES sundials_${_c}
                 HINTS "${SUNDIALS_ROOT}/lib" "${SUNDIALS_ROOT}/lib64" "${_sroot}/lib" "${_sroot}/lib64")
    if(SUNDIALS_LIB_${_c})
      list(APPEND SUNDIALS_TARGETS ${SUNDIALS_LIB_${_c}})
    elseif(NOT _c STREQUAL "core")   # sundials_core exists from SUNDIALS 7
      message(FATAL_ERROR "SUNDIALS library sundials_${_c} not found (SUNDIALS_ROOT='${SUNDIALS_ROOT}')")
    endif()
  endforeach()
  message(STATUS "SUNDIALS ${_sv}.x (library search): ${SUNDIALS_INCLUDE_DIR}")
endif()

# SuiteSparse
#
# SuiteSparse >= 7 installs one CMake package per library (KLUConfig.cmake, UMFPACKConfig.cmake, ...): point
# CMAKE_PREFIX_PATH at its prefix if needed.  Older releases are found library by library.
set(SUITESPARSE_ROOT
    "$ENV{SUITESPARSE_ROOT}"
    CACHE PATH "Optional SuiteSparse installation prefix")
set(SUITESPARSE_TARGETS "")
set(SUITESPARSE_INCLUDE_DIRS "")

# cronos_find_suitesparse(<Package> <library> <header> <dependent libraries...>)
function(cronos_find_suitesparse pkg lib header)
  find_package(${pkg} CONFIG QUIET HINTS "${SUITESPARSE_ROOT}")
  if(TARGET SuiteSparse::${pkg})
    set(_t SuiteSparse::${pkg})
  elseif(TARGET SuiteSparse::${pkg}_static)
    set(_t SuiteSparse::${pkg}_static)
  endif()
  if(_t)
    message(STATUS "  SuiteSparse ${pkg}: CMake package (${_t})")
    set(SUITESPARSE_TARGETS ${SUITESPARSE_TARGETS} ${_t} PARENT_SCOPE)
    return()
  endif()
  find_path(_inc_${pkg} ${header} HINTS "${SUITESPARSE_ROOT}/include" PATH_SUFFIXES suitesparse)
  set(_libs "")
  foreach(_l IN ITEMS ${lib} ${ARGN})
    find_library(_lib_${pkg}_${_l} NAMES ${_l} HINTS "${SUITESPARSE_ROOT}/lib" "${SUITESPARSE_ROOT}/lib64")
    if(NOT _lib_${pkg}_${_l})
      message(FATAL_ERROR "SuiteSparse library '${_l}' (for ${pkg}) not found -- set SUITESPARSE_ROOT, or turn the "
                          "corresponding CRONOS_WITH_* option off")
    endif()
    list(APPEND _libs ${_lib_${pkg}_${_l}})
  endforeach()
  if(NOT _inc_${pkg})
    message(FATAL_ERROR "SuiteSparse header '${header}' (for ${pkg}) not found -- set SUITESPARSE_ROOT")
  endif()
  message(STATUS "  SuiteSparse ${pkg}: library search (${_inc_${pkg}})")
  set(SUITESPARSE_TARGETS ${SUITESPARSE_TARGETS} ${_libs} PARENT_SCOPE)
  set(SUITESPARSE_INCLUDE_DIRS ${SUITESPARSE_INCLUDE_DIRS} ${_inc_${pkg}} PARENT_SCOPE)
endfunction()

if(CRONOS_WITH_KLU OR CRONOS_WITH_UMFPACK OR CRONOS_WITH_SPQR)
  message(STATUS "SuiteSparse:")
endif()
if(CRONOS_WITH_KLU)
  cronos_find_suitesparse(KLU klu klu.h amd colamd btf suitesparseconfig)
endif()
if(CRONOS_WITH_UMFPACK)
  cronos_find_suitesparse(UMFPACK umfpack umfpack.h amd cholmod suitesparseconfig)
endif()
if(CRONOS_WITH_SPQR)
  message(STATUS "  NOTE: CRONOS_WITH_SPQR=ON -- SPQR is GPL-2 licensed; binaries linking it are subject to the GPL")
  cronos_find_suitesparse(SPQR spqr SuiteSparseQR.hpp cholmod suitesparseconfig)
endif()

# SuperLU (through arma::spsolve: ARMA_USE_SUPERLU)
set(SUPERLU_ROOT
    "$ENV{SUPERLU_ROOT}"
    CACHE PATH "Optional SuperLU installation prefix")
if(CRONOS_WITH_SUPERLU)
  find_path(SUPERLU_INCLUDE_DIR slu_ddefs.h HINTS "${SUPERLU_ROOT}/include" PATH_SUFFIXES superlu SuperLU)
  find_library(SUPERLU_LIBRARY NAMES superlu HINTS "${SUPERLU_ROOT}/lib" "${SUPERLU_ROOT}/lib64")
  if(NOT SUPERLU_INCLUDE_DIR OR NOT SUPERLU_LIBRARY)
    message(FATAL_ERROR "SuperLU not found -- set SUPERLU_ROOT, or pass -DCRONOS_WITH_SUPERLU=OFF (OCFESLV then "
                        "falls back from SOLVE_SUPERLU)")
  endif()
  message(STATUS "SuperLU: ${SUPERLU_LIBRARY}")
endif()

# Eigen (header-only)
if(CRONOS_WITH_EIGEN)
  find_package(Eigen3 3.3 CONFIG REQUIRED)
  message(STATUS "Eigen: ${Eigen3_VERSION}")
endif()

# CRONOS Library (Interface)
#
# The Python module is named 'cronos', so the header-only library target is cronos_headers (alias
# cronos::headers): link it to build a driver against CRONOS.
add_library(cronos_headers INTERFACE)
add_library(cronos::headers ALIAS cronos_headers)

target_include_directories(
  cronos_headers
  INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
            $<INSTALL_INTERFACE:include>
            ${MCPP_INCLUDE_DIRS}
            ${Boost_INCLUDE_DIRS}
            ${ARMADILLO_INCLUDE_DIRS}
            ${INTERVAL_INCLUDE_DIRS}
            ${SUITESPARSE_INCLUDE_DIRS})
if(SUNDIALS_INCLUDE_DIR)
  target_include_directories(cronos_headers INTERFACE ${SUNDIALS_INCLUDE_DIR})
endif()
if(CRONOS_WITH_SUPERLU)
  target_include_directories(cronos_headers INTERFACE ${SUPERLU_INCLUDE_DIR})
endif()

# Compile definitions: MC++'s, then CRONOS's backends
target_compile_definitions(
  cronos_headers INTERFACE MC__USE_THREAD MC__USE_ARMADILLO ${INTERVAL_DEFINITIONS})
if(MC__USE_FADBAD)
  target_compile_definitions(cronos_headers INTERFACE MC__USE_FADBAD)
endif()
if(ENABLE_HSL)
  target_compile_definitions(cronos_headers INTERFACE MC__USE_HSL)
endif()
if(CRONOS_WITH_KLU)
  target_compile_definitions(cronos_headers INTERFACE CRONOS__WITH_KLU)
endif()
if(CRONOS_WITH_UMFPACK)
  target_compile_definitions(cronos_headers INTERFACE CRONOS__WITH_UMFPACK)
endif()
if(CRONOS_WITH_SPQR)
  target_compile_definitions(cronos_headers INTERFACE CRONOS__WITH_SPQR)
endif()
if(CRONOS_WITH_EIGEN)
  target_compile_definitions(cronos_headers INTERFACE CRONOS__WITH_EIGEN)
endif()
if(CRONOS_WITH_SUPERLU)
  target_compile_definitions(cronos_headers INTERFACE ARMA_USE_SUPERLU)
endif()

# Match the classical Make warning flags and add backend-specific compile options such as -frounding-math only
# when required, e.g. with FILIB.
if(MSVC)
  target_compile_options(cronos_headers INTERFACE ${INTERVAL_COMPILE_OPTIONS})
else()
  target_compile_options(
    cronos_headers
    INTERFACE -Wall
              -Wno-misleading-indentation
              -Wno-unknown-pragmas
              -Wno-parentheses
              -Wno-unused-result
              -Wno-unused-function
              -Wno-comment
              ${INTERVAL_COMPILE_OPTIONS})
endif()

# Link libraries (SuperLU before LAPACK/BLAS, which it uses)
if(CRONOS_WITH_SUPERLU)
  target_link_libraries(cronos_headers INTERFACE ${SUPERLU_LIBRARY})
endif()
target_link_libraries(
  cronos_headers
  INTERFACE ${SUNDIALS_TARGETS}
            ${SUITESPARSE_TARGETS}
            ${ARMADILLO_LIBRARIES}
            ${LAPACK_LIBRARIES}
            ${BLAS_LIBRARIES}
            ${INTERVAL_LIBRARIES}
            Threads::Threads
            ${CMAKE_DL_LIBS})   # dlsym: run-time OpenMP / OpenBLAS thread control in ffmodel.hpp
if(CRONOS_WITH_EIGEN)
  target_link_libraries(cronos_headers INTERFACE Eigen3::Eigen)
endif()
if(ENABLE_HSL)
  target_link_libraries(cronos_headers INTERFACE ${HSL_LIBRARIES})
endif()

# Python Library
if(ENABLE_PYTHON)
  set(CRONOS_PY_SOURCES
      src/interface/ffdom.cpp
      src/interface/ffmodel.cpp
      src/interface/odeslv.cpp
      src/interface/ffode.cpp
      src/interface/ocfeslv.cpp
      src/interface/ffocfe.cpp
      src/interface/cronos.cpp)

  pybind11_add_module(cronos ${CRONOS_PY_SOURCES})

  if(CRONOS_VENDORED_PYBIND11)
    target_include_directories(
      cronos BEFORE PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/extern/pybind11/include")
  endif()
  # the generated option headers (gen_*_options.hpp) live next to the binders
  target_include_directories(cronos PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src/interface")
  target_link_libraries(cronos PRIVATE cronos_headers)

  # Type stub: auto-generated at build time, as for pymcpp.  pybind11-stubgen IMPORTS the freshly-built module,
  # which imports pymcpp (FFVar, FFGraph, ... are pymcpp's types), so both pybind11-stubgen and pymcpp must be
  # importable by the build interpreter; PYMCPP_DIR is the directory holding the pymcpp module, if it is not on the
  # interpreter's path.  This step cannot run when cross-compiling: pass -DCRONOS_STUBS=OFF then.
  option(CRONOS_STUBS "Generate and install the .pyi type stub" ON)
  set(PYMCPP_DIR
      ""
      CACHE PATH "Directory holding the pymcpp module (for stub generation), if not on the Python path")

  # PYTHONPATH entries are joined with ";" on Windows and ":" elsewhere.
  if(WIN32)
    set(_stub_pathsep "$<SEMICOLON>")
    set(_probe_pathsep ";")
  else()
    set(_stub_pathsep ":")
    set(_probe_pathsep ":")
  endif()

  # CRONOS_STUBS records what the user ASKED for; probe the build interpreter at configure time and keep the
  # verdict apart, so a missing pybind11-stubgen or pymcpp does not fail the build after cronos has linked.
  set(_cronos_stubs_ok ${CRONOS_STUBS})
  if(CRONOS_STUBS)
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E env "PYTHONPATH=${PYMCPP_DIR}${_probe_pathsep}$ENV{PYTHONPATH}"
              "${Python3_EXECUTABLE}" -c "import pybind11_stubgen, pymcpp"
      RESULT_VARIABLE _cronos_stubgen_rc
      OUTPUT_QUIET ERROR_QUIET)
    if(NOT _cronos_stubgen_rc EQUAL 0)
      set(_cronos_stubs_ok OFF)
      message(WARNING
        "pybind11-stubgen and pymcpp are not both importable by ${Python3_EXECUTABLE} -- the cronos.pyi type stub "
        "will not be generated or installed.  The extension module itself is unaffected.  To enable it:\n"
        "    pip install pybind11-stubgen\n"
        "    and set -DPYMCPP_DIR=<directory holding the pymcpp module> if pymcpp is not on the Python path.\n"
        "Pass -DCRONOS_STUBS=OFF to silence this warning.")
    endif()
  endif()

  if(_cronos_stubs_ok)
    set(_stub_out ${CMAKE_CURRENT_BINARY_DIR}/stubs)
    add_custom_command(
      TARGET cronos
      POST_BUILD
      COMMAND
        ${CMAKE_COMMAND} -E env
        "PYTHONPATH=$<TARGET_FILE_DIR:cronos>${_stub_pathsep}${PYMCPP_DIR}${_stub_pathsep}$ENV{PYTHONPATH}"
        ${Python3_EXECUTABLE} -m pybind11_stubgen cronos -o ${_stub_out}
      BYPRODUCTS ${_stub_out}/cronos.pyi
      COMMENT "Generating cronos.pyi type stub"
      VERBATIM)
  endif()
endif()

# `make check` (or `cmake --build <dir> --target check`): runs the five Python check scripts against the freshly
# built cronos module, in order, stopping at the first failure (a script exits non-zero when one of its checks
# fails).  The scripts find the modules through CRONOS_PYPATH: the cronos build directory, then PYMCPP_DIR.
if(ENABLE_PYTHON)
  set(CRONOS_PYCHECK_DIR
      "${CMAKE_CURRENT_SOURCE_DIR}/src/interface"
      CACHE PATH "Directory holding the Python check scripts (check_*.py)")
  set(CRONOS_PYCHECKS
      check_ffmodel.py
      check_odeslv.py
      check_ffodeslv.py
      check_ocfeslv.py
      check_ffocfe.py)
  set(_check_cmds "")
  set(_check_missing "")
  foreach(_chk IN LISTS CRONOS_PYCHECKS)
    if(EXISTS "${CRONOS_PYCHECK_DIR}/${_chk}")
      list(APPEND _check_cmds
           COMMAND ${CMAKE_COMMAND} -E echo "== ${_chk}"
           COMMAND ${CMAKE_COMMAND} -E env
                   "CRONOS_PYPATH=$<TARGET_FILE_DIR:cronos>${_stub_pathsep}${PYMCPP_DIR}"
                   "MPLBACKEND=Agg"
                   ${Python3_EXECUTABLE} "${CRONOS_PYCHECK_DIR}/${_chk}")
    else()
      list(APPEND _check_missing ${_chk})
    endif()
  endforeach()
  if(_check_missing)
    message(WARNING "make check: not found in CRONOS_PYCHECK_DIR (${CRONOS_PYCHECK_DIR}): ${_check_missing}")
  endif()
  if(_check_cmds)
    add_custom_target(
      check
      ${_check_cmds}
      DEPENDS cronos
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMENT "Running the cronos Python checks"
      USES_TERMINAL VERBATIM)
  endif()
endif()

# Installation
#
# Direct CMake installs default to: ${CMAKE_INSTALL_PREFIX}/include for the CRONOS headers,
# ${CMAKE_INSTALL_PREFIX}/lib for the compiled cronos extension module, ${CMAKE_INSTALL_PREFIX}/notebook for the
# notebooks and Python scripts.  For wheel builds with scikit-build-core, override CRONOS_INSTALL_DIR with
# -DCRONOS_INSTALL_DIR=. and select only the python_modules component.
set(CRONOS_INSTALL_DIR
    "lib"
    CACHE STRING "Installation directory for the cronos module")
set(CRONOS_NOTEBOOK_INSTALL_DIR
    "notebook"
    CACHE STRING "Installation directory for notebooks and Python scripts")

message(STATUS "--------------------------------------------------------")
message(STATUS "Backends: KLU=${CRONOS_WITH_KLU} UMFPACK=${CRONOS_WITH_UMFPACK} SPQR=${CRONOS_WITH_SPQR} "
               "Eigen=${CRONOS_WITH_EIGEN} SuperLU=${CRONOS_WITH_SUPERLU} HSL=${ENABLE_HSL} FADBAD=${MC__USE_FADBAD}")
message(
  STATUS "Installation Prefix (CMAKE_INSTALL_PREFIX): ${CMAKE_INSTALL_PREFIX}")
message(STATUS "Headers will be installed to: ${CMAKE_INSTALL_PREFIX}/include")
if(ENABLE_PYTHON)
  message(
    STATUS
      "cronos module will be installed to: ${CMAKE_INSTALL_PREFIX}/${CRONOS_INSTALL_DIR}")
endif()
message(
  STATUS
    "Notebooks/scripts will be installed to: ${CMAKE_INSTALL_PREFIX}/${CRONOS_NOTEBOOK_INSTALL_DIR}"
)
message(STATUS "--------------------------------------------------------")

# Install headers
set(CRONOS_HEADERS
    base_cvodes.hpp
    base_sundials.hpp
    ffmodel.hpp
    ffocfe.hpp
    ffode.hpp
    ocfeslv.hpp
    ocplan.hpp
    odeslv_base.hpp
    odeslv_cvodes.hpp
    odeslvs_base.hpp
    odeslvs_cvodes.hpp)

set(CRONOS_HEADERS_FULL "")
foreach(HEADER ${CRONOS_HEADERS})
  list(APPEND CRONOS_HEADERS_FULL "src/${HEADER}")
endforeach()

install(
  FILES ${CRONOS_HEADERS_FULL}
  DESTINATION include
  COMPONENT headers)

if(ENABLE_PYTHON)
  install(
    TARGETS cronos
    DESTINATION ${CRONOS_INSTALL_DIR}
    COMPONENT python_modules)
  # Tests the configure-time verdict, not the option, so a missing pybind11-stubgen does not leave install()
  # looking for a stub that was never generated.
  if(_cronos_stubs_ok)
    install(
      FILES ${CMAKE_CURRENT_BINARY_DIR}/stubs/cronos.pyi
      DESTINATION ${CRONOS_INSTALL_DIR}
      COMPONENT python_modules)
  endif()
endif()

# Install notebooks and Python scripts, if the notebook/ directory exists.
install(
  DIRECTORY notebook/
  DESTINATION ${CRONOS_NOTEBOOK_INSTALL_DIR}
  COMPONENT notebooks
  OPTIONAL FILES_MATCHING
  PATTERN "*.ipynb"
  PATTERN "*.py"
  PATTERN "*.md")

# Examples: the test drivers, one executable each (target X from X.cpp, as the Make setup), registered with CTest
# (a driver exits non-zero on failure).  Labels: the directory (ODESLV, OCFESLV) and "slow" for the drivers taking
# 30 s or more (measured on the 3 Oct 2026 sweep) -- `ctest -LE slow` runs the quick set.
set(CRONOS_SLOW_TESTS
    OCFE_PDE5 OCFE_PDE5b OCFE_PSA7 OCFE_PSA10 OCFE_MBC3 OCFE_MBC4 OCFE_MBC5 OCFE_MBC5b OCFE_MBC5c
    OCFE_MBC5fig OCFE_MBC5fig_ebal OCFE_MBC6 OCFE_MBC6fig OCFE_MBC6fig_ebal
    CACHE STRING "Test drivers labelled 'slow' (>= 30 s)")
set(CRONOS_TEST_TIMEOUT 300 CACHE STRING "CTest timeout (s) for a driver")
set(CRONOS_SLOW_TEST_TIMEOUT 1800 CACHE STRING "CTest timeout (s) for a 'slow' driver")
if(ENABLE_EXAMPLES)
  enable_testing()
  # the ODESLV drivers' directory has been named both ways
  foreach(_dir ODESLV ODESLVS OCFESLV)
    if(NOT IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/test/${_dir}")
      continue()
    endif()
    file(GLOB _drivers CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/test/${_dir}/*.cpp")
    foreach(SOURCE_FILE IN LISTS _drivers)
      get_filename_component(EXEC_NAME "${SOURCE_FILE}" NAME_WE)
      add_executable("${EXEC_NAME}" "${SOURCE_FILE}")
      target_link_libraries("${EXEC_NAME}" PRIVATE cronos_headers)
      # helpers shared by the drivers of a directory (e.g. test_deriv_utils.hpp)
      target_include_directories("${EXEC_NAME}" PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/test/${_dir}")
      if(_dir STREQUAL "OCFESLV")
        target_compile_definitions("${EXEC_NAME}" PRIVATE OCFE_OCFESLV_HEADER="ocfeslv.hpp")
      endif()
      add_test(NAME "${EXEC_NAME}" COMMAND "${EXEC_NAME}"
               WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
      if("${EXEC_NAME}" IN_LIST CRONOS_SLOW_TESTS)
        set_tests_properties("${EXEC_NAME}" PROPERTIES LABELS "${_dir};slow" TIMEOUT ${CRONOS_SLOW_TEST_TIMEOUT})
      else()
        set_tests_properties("${EXEC_NAME}" PROPERTIES LABELS "${_dir}" TIMEOUT ${CRONOS_TEST_TIMEOUT})
      endif()
    endforeach()
  endforeach()
endif()
