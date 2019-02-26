# Third-Party Libraries -->> CHANGE ACCORDINGLY <<--
#PROF = -pg
#OPTI = -Ofast
#OPTI = -O2

PATH_MC = $(HOME)/Programs/Devel/MC++/MC++_2.1
LIB_MC = -llapack -lblas -lmc13 -lmc21 -lmc33 -lgfortran
INC_MC = -I$(PATH_MC)/src/mc -I$(PATH_MC)/src/3rdparty/fadbad++ -I$(PATH_MC)/src/3rdparty/cpplapack-2015.05.11-1/include
FLAGS_MC = -Wno-misleading-indentation -Wno-unknown-pragmas -DMC__USE_HSL

PATH_PROFIL = /home/bchachua/Programs/ThirdParty/Profil-2.0.8
LIB_PROFIL  = -L$(PATH_PROFIL)/lib -lProfilPackages -lProfil -lBias -llr
INC_PROFIL  = -I$(PATH_PROFIL)/include

PATH_FILIB  = /opt/filib++
LIB_FILIB   = -L$(PATH_FILIB)/lib -lprim
INC_FILIB   = -I$(PATH_FILIB)/include/ -I$(PATH_FILIB)/include/interval
FLAGS_FILIB = -frounding-math -ffloat-store

PATH_IPOPT   = /opt/CoinIpopt
PATH_COINHSL = /opt/coinhsl
LIB_IPOPT = -L$(PATH_IPOPT)/lib -L$(PATH_COINHSL)/lib -lipopt -lcoinhsl -llapack -lblas -lgfortran -lm -lquadmath -lblas -lgfortran -lm -lquadmath -lm -ldl -lcoinmumps -lblas -lgfortran -lm -lquadmath -lgfortran -lm -lquadmath -lcoinmetis
INC_IPOPT = -I$(PATH_IPOPT)/include

PATH_CPLEX   = /opt/ibm/ILOG/CPLEX_Studio1263/cplex
PATH_CONCERT = /opt/ibm/ILOG/CPLEX_Studio1263/concert
LIB_MIP      = -L$(PATH_CPLEX)/lib/x86-64_linux/static_pic -lilocplex -lcplex \
               -L$(PATH_CONCERT)/lib/x86-64_linux/static_pic -lconcert \
               -lm -pthread
INC_MIP      = -I$(PATH_CPLEX)/include -I$(PATH_CONCERT)/include
FLAGS_MIP    = -DMC__USE_CPLEX -m64 -fPIC -fexceptions -DIL_STD -Wno-ignored-attributes

#PATH_GUROBI = $(GUROBI_HOME)
#LIB_MIP     = -L$(PATH_GUROBI)/lib -lgurobi_g++5.2 -lgurobi75 -pthread
#INC_MIP     = -I$(PATH_GUROBI)/include
#FLAGS_MIP   = -DMC__USE_GUROBI

FLAG_ALL = -g $(OPTI) -Wall -std=c++11 $(FLAGS_FILIB) $(FLAGS_MC) $(FLAGS_MIP)
LIB_ALL  = $(LIB_MC) $(LIB_PROFIL) $(LIB_FILIB) $(LIB_MIP) $(LIB_IPOPT)
INC_ALL  = $(INC_MC) $(INC_PROFIL) $(INC_FILIB) $(INC_MIP) $(INC_IPOPT)

#####

all: test0 test1 test2 test3 test4 test5 test6 test7 test8 test9

test% : test%.o
	g++ $(PROF) -o $@ $< $(LIB_ALL)

test%.o: test%.cpp
	g++ $(PROF) -c $(FLAG_ALL) -I../../src $(INC_ALL) -o $@ $<

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
	rm -f test? test?? *.o *~ *.eps *.lp *.out *.log
