# THIRD-PARTY LIBRARIES <<-- CHANGE AS APPROPRIATE -->>

PATH_CRONOS_MK := $(dir $(lastword $(MAKEFILE_LIST)))

include $(abspath $(PATH_CRONOS_MK)../../mcpp/src/makeoptions.mk)
PATH_MC = $(abspath $(PATH_CRONOS_MK)../../mcpp)

PATH_CRONOS = $(abspath ../)

PATH_SUNDIALS = $(SUNDIALS_HOME)
LIB_SUNDIALS = -L$(PATH_SUNDIALS)/lib -lsundials_sunlinsolklu -lsundials_cvodes -lsundials_nvecserial -lsundials_core -llapack -lblas
INC_SUNDIALS = -I/usr/include/suitesparse -I$(PATH_SUNDIALS)/include
FLAG_SUNDIALS = -DCRONOS__WITH_KLU

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

FLAG_CRONOS = -fPIC $(FLAG_MC) $(FLAG_SUNDIALS)
LIB_CRONOS  = $(LIB_MC) $(LIB_SUNDIALS)
INC_CRONOS  = -I$(PATH_CRONOS)/src $(INC_MC) $(INC_SUNDIALS)

