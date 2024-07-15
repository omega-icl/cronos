# THIRD-PARTY LIBRARIES <<-- CHANGE AS APPROPRIATE -->>

PATH_MC = $(shell cd $(HOME)/Programs/bitbucket/mcpp30 ; pwd)
include $(PATH_MC)/src/makeoptions.mk

PATH_CRONOS    = $(shell cd $(HOME)/Programs/bitbucket/cronos; pwd)

PATH_SUNDIALS = /opt/sundials-7.0.0
LIB_SUNDIALS = -L$(PATH_SUNDIALS)/lib -lsundials_cvodes -lsundials_nvecserial -lsundials_core -llapack -lblas
INC_SUNDIALS = -I$(PATH_SUNDIALS)/include
FLAG_SUNDIALS =

#INC_PYBIND11 = -I/usr/include/python3.10 -I$(PATH_EXTERN)/pybind11/include

# COMPILATION <<-- CHANGE AS APPROPRIATE -->>

PROF = #-pg
OPTIM = -O2
DEBUG = #-g
WARN  = -Wall -Wno-misleading-indentation -Wno-unknown-pragmas -Wno-parentheses -Wno-unused-result
CPP17 = -std=c++17
CC    = gcc-13
CPP   = g++-13
# CPP   = icpc

# <<-- NO CHANGE BEYOND THIS POINT -->>

FLAG_CPP  = $(DEBUG) $(OPTIM) $(CPP17) $(WARN) $(PROF)
LINK      = $(CPP)
FLAG_LINK = $(PROF)

FLAG_CRONOS = -fPIC $(FLAG_MC) $(FLAG_SUNDIALS) 
LIB_CRONOS  = $(LIB_MC) $(LIB_SUNDIALS)
INC_CRONOS  = -I$(PATH_CRONOS)/src $(INC_MC) $(INC_SUNDIALS)

