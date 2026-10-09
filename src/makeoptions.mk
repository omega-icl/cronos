# THIRD-PARTY LIBRARIES <<-- CHANGE AS APPROPRIATE -->>

PATH_CRONOS_MK := $(dir $(lastword $(MAKEFILE_LIST)))

#include $(abspath $(PATH_CRONOS_MK)../../mcpp/src/makeoptions.mk)
#PATH_MC = $(abspath $(PATH_CRONOS_MK)../../mcpp)
PATH_MC = /home/bchachua/Programs/github/mcpp
include $(abspath $(PATH_MC)/src/makeoptions.mk)

PATH_CRONOS = $(abspath $(PATH_CRONOS_MK)../)

PATH_SUITESPARSE = 
LIB_SUITESPARSE = -lspqr -lcholmod -lumfpack -lsuitesparseconfig -lamd -lcolamd -lccolamd -lcamd -lmetis
INC_SUITESPARSE = -I/usr/include/suitesparse
FLAG_SUITESPARSE = -DCRONOS__WITH_KLU -DCRONOS__WITH_SPQR -DCRONOS__WITH_UMFPACK

PATH_EIGEN = 
LIB_EIGEN = 
INC_EIGEN = -I/usr/include/eigen3
FLAG_EIGEN = -DCRONOS__WITH_EIGEN

PATH_SUPERLU = 
LIB_SUPERLU = -lsuperlu
INC_SUPERLU = -I/usr/include/superlu
FLAG_SUPERLU =

PATH_SUNDIALS = $(SUNDIALS_HOME)
LIB_SUNDIALS = -L$(PATH_SUNDIALS)/lib -lsundials_sunlinsolklu -lsundials_cvodes -lsundials_nvecserial -lsundials_core -llapack -lblas
INC_SUNDIALS = -I$(PATH_SUNDIALS)/include
FLAG_SUNDIALS = 


# COMPILATION <<-- CHANGE AS APPROPRIATE -->>

PROF = #-pg
OPTIM = -O2
DEBUG = #-g
WARN  = -Wall -Wno-misleading-indentation -Wno-unknown-pragmas -Wno-parentheses -Wno-unused-result
CPP17 = -std=c++17
CC    = gcc
CPP   = g++
# CPP   = icpc

# <<-- NO CHANGE BEYOND THIS POINT -->>

FLAG_CPP  = $(DEBUG) $(OPTIM) $(CPP17) $(WARN) $(PROF)
LINK      = $(CPP)
FLAG_LINK = $(PROF)

FLAG_CRONOS = -fPIC $(FLAG_MC) $(FLAG_SUITESPARSE) $(FLAG_EIGEN) $(FLAG_SUNDIALS)
LIB_CRONOS  = $(LIB_MC) $(LIB_SUITESPARSE) $(LIB_EIGEN) $(LIB_SUNDIALS)
INC_CRONOS  = -I$(PATH_CRONOS)/src $(INC_MC) $(INC_SUITESPARSE) $(INC_EIGEN) $(INC_SUNDIALS)

