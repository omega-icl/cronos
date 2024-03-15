# THIRD-PARTY LIBRARIES <<-- CHANGE AS APPROPRIATE -->>

PATH_MC = $(shell cd $(HOME)/Programs/bitbucket/mcpp30 ; pwd)
include $(PATH_MC)/src/makeoptions.mk

PATH_SUNDIALS = /opt/sundials-7.0.0
LIB_SUNDIALS = -L$(PATH_SUNDIALS)/lib -lsundials_cvodes -lsundials_nvecserial -lsundials_core -llapack -lblas
INC_SUNDIALS = -I$(PATH_SUNDIALS)/include
FLAG_SUNDIALS =

FLAG_DEP = -fPIC $(FLAG_MC) $(FLAG_SUNDIALS) 
LIB_DEP  = $(LIB_MC) $(LIB_SUNDIALS)
INC_DEP  = $(INC_MC) $(INC_SUNDIALS)

# COMPILATION <<-- CHANGE AS APPROPRIATE -->>

DEBUG = -g
#PROF = -pg
OPTIM = #-O2 #-Ofast
WARN  = -Wall -Wno-misleading-indentation -Wno-unknown-pragmas -Wno-unused-result
CPP17 = -std=c++17

CC  = gcc-12
CPP = g++-12
# CPP = icpc
FLAG_CPP = $(DEBUG) $(PROF) $(OPTIM) $(CPP17) $(WARN) $(FLAG_DEP) 

LINK = $(CPP)
FLAG_LINK = 

#LDFLAGS = -ldl -Wl,-rpath,\$$ORIGIN -Wl,-rpath,$(PATH_GAMS)

