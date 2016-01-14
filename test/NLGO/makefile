# Third-Party Libraries -->> CHANGE ACCORDINGLY <<--
PATH_MC = /home/bchachua/Programs/Devel/MC++/MC++_2.0
LIB_MC = -llapack -lblas
INC_MC = -I$(PATH_MC)/include -I$(PATH_MC)/src/3rdparty/fadbad++ -I$(PATH_MC)/src/3rdparty/cpplapack/include

LIB_HSL = -lmc13 -lmc21
INC_HSL = 

PATH_PROFIL = /home/bchachua/Programs/ThirdParty/Profil-2.0.8
LIB_PROFIL  = -L$(PATH_PROFIL)/lib -lProfilPackages -lProfil -lBias -llr
INC_PROFIL  = -I$(PATH_PROFIL)/include

PATH_FILIB  = /opt/filib++
LIB_FILIB   = -L$(PATH_FILIB)/lib -lprim
INC_FILIB   = -I$(PATH_FILIB)/include/ -I$(PATH_FILIB)/include/interval
FLAGS_FILIB = -frounding-math -ffloat-store

PATH_IPOPT  = /opt/Ipopt-3.10.2
LIB_IPOPT   = -L$(PATH_IPOPT)/lib -lipopt -lcoinhsl -llapack -lblas -ldl
INC_IPOPT   = -I$(PATH_IPOPT)/include

PATH_CPLEX   = /opt/ibm/ILOG/CPLEX_Studio1251/cplex
PATH_CONCERT = /opt/ibm/ILOG/CPLEX_Studio1251/concert
LIB_CPLEX    = -L$(PATH_CPLEX)/lib/x86-64_sles10_4.1/static_pic -lilocplex -lcplex \
               -L$(PATH_CONCERT)/lib/x86-64_sles10_4.1/static_pic -lconcert \
	       -lm -pthread
INC_CPLEX    = -I$(PATH_CPLEX)/include -I$(PATH_CONCERT)/include
FLAGS_CPLEX  = -m64 -fPIC -fexceptions -DIL_STD
#FLAGS_CPLEX  = -m64 -fPIC -fexceptions -DNDEBUG -DIL_STD

PATH_GUROBI = $(GUROBI_HOME)
LIB_GUROBI  = -L$(PATH_GUROBI)/lib -lgurobi_c++ -lgurobi56 -pthread
INC_GUROBI  = -I$(PATH_GUROBI)/include

#LIB_SDPA   = -lsdpa -ldmumps_seq
#INC_SDPA   = 

#PROF = -pg

FLAG_ALL = $(PROF) -g -O2 -Wall -std=c++11 $(FLAG_FILIB)
LIB_ALL  = $(LIB_MC) $(LIB_PROFIL) $(LIB_FILIB) $(LIB_HSL) $(LIB_LAPACK) $(LIB_GUROBI) $(LIB_CPLEX) $(LIB_IPOPT)
INC_ALL  = $(INC_MC) $(INC_PROFIL) $(INC_FILIB) $(INC_HSL) $(INC_LAPACK) $(INC_GUROBI) $(INC_CPLEX) $(INC_IPOPT)

#####

all: test0 test1 test2 test3 test4 test4b test5 test5b test6 test7 test8 test9

test% : test%.o
	g++ -o $@ $< $(LIB_ALL)

test%.o: test%.cpp
	g++ -c $(FLAG_ALL) -I../../src $(INC_ALL) -o $@ $<

#####

graph : gendot gensvg

gendot :
	dot -Tsvg -O $(shell ls *.dot)

gensvg :
	@for FILE in $(shell ls *.svg); do \
	  echo $$FILE; \
	  (display $$FILE &); \
	done

#####

clean :
	rm -f test? test?b *.o *~ *.eps *.lp

