# THIRD-PARTY LIBRARIES <<-- CHANGE AS APPROPRIATE -->>

PATH_MC = $(HOME)/Programs/Devel/MC++/MC++_2.1
LIB_MC = -llapack -lblas -lmc13 -lmc21 -lmc33 -lgfortran
INC_MC = -I$(PATH_MC)/src/mc -I$(PATH_MC)/src/3rdparty/fadbad++ -I$(PATH_MC)/src/3rdparty/cpplapack-2015.05.11-1/include
FLAG_MC = -Wno-misleading-indentation -Wno-unknown-pragmas -DMC__USE_HSL

PATH_PROFIL = /home/bchachua/Programs/ThirdParty/Profil-2.0.8
LIB_PROFIL  = -L$(PATH_PROFIL)/lib -lProfilPackages -lProfil -lBias -llr
INC_PROFIL  = -I$(PATH_PROFIL)/include

PATH_FILIB  = /opt/filib++
LIB_FILIB   = -L$(PATH_FILIB)/lib -lprim
INC_FILIB   = -I$(PATH_FILIB)/include -I$(PATH_FILIB)/include/interval
FLAG_FILIB = -frounding-math -ffloat-store

PATH_SOBOL  = $(HOME)/Programs/Devel/CRONOS/CRONOS_1.1/src/3rdparty
LIB_SOBOL   = -L$(PATH_SOBOL) -lsobol
INC_SOBOL   = -I$(PATH_SOBOL)
FLAG_SOBOL  = -DMC__USE_SOBOL

PATH_IPOPT   = /opt/CoinIpopt
PATH_COINHSL = /opt/coinhsl
LIB_IPOPT = -L$(PATH_IPOPT)/lib -L$(PATH_COINHSL)/lib -lipopt -lcoinhsl -llapack -lblas -lgfortran -lm -lquadmath -lblas -lgfortran -lm -lquadmath -lm -ldl -lcoinmumps -lblas -lgfortran -lm -lquadmath -lgfortran -lm -lquadmath -lcoinmetis
INC_IPOPT = -I$(PATH_IPOPT)/include

PATH_SUNDIALS = /opt/sundials-2.7.0
LIB_SUNDIALS = -L$(PATH_SUNDIALS)/lib -lsundials_cvodes -lsundials_nvecserial -llapack -lblas
INC_SUNDIALS = -I$(PATH_SUNDIALS)/include -I$(PATH_SUNDIALS)/include/sundials -I$(PATH_SUNDIALS)/include/cvodes -I$(PATH_SUNDIALS)/include/nvector

PATH_CPLEX   = /opt/ibm/ILOG/CPLEX_Studio1263/cplex
PATH_CONCERT = /opt/ibm/ILOG/CPLEX_Studio1263/concert
LIB_MIP      = -L$(PATH_CPLEX)/lib/x86-64_linux/static_pic -lilocplex -lcplex \
               -L$(PATH_CONCERT)/lib/x86-64_linux/static_pic -lconcert \
               -lm -pthread
INC_MIP      = -I$(PATH_CPLEX)/include -I$(PATH_CONCERT)/include
FLAG_MIP    = -DMC__USE_CPLEX -m64 -fPIC -fexceptions -DIL_STD -Wno-ignored-attributes

#PATH_GUROBI = $(GUROBI_HOME)
#LIB_MIP     = -L$(PATH_GUROBI)/lib -lgurobi_g++5.2 -lgurobi75 -pthread
#INC_MIP     = -I$(PATH_GUROBI)/include
#FLAG_MIP   = -DMC__USE_GUROBI

FLAG_DEP = $(FLAG_MC) $(FLAG_FILIB) $(FLAG_MIP) $(FLAG_SOBOL) 
LIB_DEP  = $(LIB_MC) $(LIB_PROFIL) $(LIB_FILIB) $(LIB_MIP) $(LIB_IPOPT) $(LIB_SUNDIALS) $(LIB_SOBOL)
INC_DEP  = $(INC_MC) $(INC_PROFIL) $(INC_FILIB) $(INC_MIP) $(INC_IPOPT) $(INC_SUNDIALS) $(INC_SOBOL)

# COMPILATION <<-- CHANGE AS APPROPRIATE -->>

DEBUG = -g
#PROF = -pg
OPTIM = -Ofast
WARN  = -Wall
CPP17 = -std=c++1z

CPP = g++-7
# CPP = icpc
FLAG_CPP = $(DEBUG) $(PROF) $(OPTIM) $(CPP17) $(WARN) $(FLAG_DEP) 

LINK = $(CPP)
FLAG_LINK = 

