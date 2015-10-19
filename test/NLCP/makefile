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

#####

all: test1 test2 test3 test4 test5 test5b test6 test7 test8 test9

test1 : test1.o
	g++ $(PROF) -o test1 test1.o $(LIB_MC) $(LIB_PROFIL) $(LIB_FILIB) $(LIB_HSL) $(LIB_LAPACK) $(LIB_GUROBI) $(LIB_IPOPT)

test1.o: test1.cpp
	g++ -c $(PROF) -g -Wall -std=c++11 $(FLAGS_FILIB) -I../../src $(INC_MC) $(INC_PROFIL) $(INC_FILIB) $(INC_LAPACK) $(INC_GUROBI) $(INC_IPOPT) -o test1.o test1.cpp

test2 : test2.o
	g++ $(PROF) -o test2 test2.o $(LIB_MC) $(LIB_PROFIL) $(LIB_FILIB) $(LIB_HSL) $(LIB_LAPACK) $(LIB_GUROBI) $(LIB_IPOPT)

test2.o: test2.cpp
	g++ -c $(PROF) -g -Wall -std=c++11 $(FLAGS_FILIB) -I../../src $(INC_MC) $(INC_PROFIL) $(INC_FILIB) $(INC_LAPACK) $(INC_GUROBI) $(INC_IPOPT) -o test2.o test2.cpp

test3 : test3.o
	g++ $(PROF) -o test3 test3.o $(LIB_MC) $(LIB_PROFIL) $(LIB_FILIB) $(LIB_HSL) $(LIB_LAPACK) $(LIB_GUROBI) $(LIB_IPOPT)

test3.o: test3.cpp
	g++ -c $(PROF) -g -Wall -std=c++11 $(FLAGS_FILIB) -I../../src $(INC_MC) $(INC_PROFIL) $(INC_FILIB) $(INC_LAPACK) $(INC_GUROBI) $(INC_IPOPT) -o test3.o test3.cpp

test4 : test4.o
	g++ $(PROF) -o test4 test4.o $(LIB_MC) $(LIB_PROFIL) $(LIB_FILIB) $(LIB_HSL) $(LIB_LAPACK) $(LIB_GUROBI) $(LIB_IPOPT)

test4.o: test4.cpp
	g++ -c $(PROF) -g -O2 -Wall -std=c++11 $(FLAGS_FILIB) -I../../src $(INC_MC) $(INC_PROFIL) $(INC_FILIB) $(INC_LAPACK) $(INC_GUROBI) $(INC_IPOPT) -o test4.o test4.cpp

test5 : test5.o
	g++ $(PROF) -o test5 test5.o $(LIB_MC) $(LIB_PROFIL) $(LIB_FILIB) $(LIB_HSL) $(LIB_LAPACK) $(LIB_GUROBI) $(LIB_IPOPT)

test5.o: test5.cpp
	g++ -c $(PROF) -g -O2 -Wall -std=c++11 $(FLAGS_FILIB) -I../../src $(INC_MC) $(INC_PROFIL) $(INC_FILIB) $(INC_LAPACK) $(INC_GUROBI) $(INC_IPOPT) -o test5.o test5.cpp

test5b : test5b.o
	g++ $(PROF) -o test5b test5b.o $(LIB_MC) $(LIB_PROFIL) $(LIB_FILIB) $(LIB_HSL) $(LIB_LAPACK) $(LIB_GUROBI) $(LIB_IPOPT)

test5b.o: test5b.cpp
	g++ -c $(PROF) -g -O2 -Wall -std=c++11 $(FLAGS_FILIB) -I../../src $(INC_MC) $(INC_PROFIL) $(INC_FILIB) $(INC_LAPACK) $(INC_GUROBI) $(INC_IPOPT) -o test5b.o test5b.cpp

test6 : test6.o
	g++ $(PROF) -o test6 test6.o $(LIB_MC) $(LIB_PROFIL) $(LIB_FILIB) $(LIB_HSL) $(LIB_LAPACK) $(LIB_GUROBI) $(LIB_IPOPT)

test6.o: test6.cpp
	g++ -c $(PROF) -g -O2 -Wall -std=c++11 $(FLAGS_FILIB) -I../../src $(INC_MC) $(INC_PROFIL) $(INC_FILIB) $(INC_LAPACK) $(INC_GUROBI) $(INC_IPOPT) -o test6.o test6.cpp

test6b : test6b.o
	g++ $(PROF) -o test6b test6b.o $(LIB_MC) $(LIB_PROFIL) $(LIB_FILIB) $(LIB_HSL) $(LIB_LAPACK) $(LIB_GUROBI) $(LIB_IPOPT)

test6b.o: test6b.cpp
	g++ -c $(PROF) -g -O2 -Wall -std=c++11 $(FLAGS_FILIB) -I../../src $(INC_MC) $(INC_PROFIL) $(INC_FILIB) $(INC_LAPACK) $(INC_GUROBI) $(INC_IPOPT) -o test6b.o test6b.cpp

test7 : test7.o
	g++ $(PROF) -o test7 test7.o $(LIB_MC) $(LIB_PROFIL) $(LIB_FILIB) $(LIB_HSL) $(LIB_LAPACK) $(LIB_GUROBI) $(LIB_CPLEX) $(LIB_IPOPT)

test7.o: test7.cpp
	g++ -c $(PROF) -g -O2 -Wall -std=c++11 $(FLAGS_FILIB) $(FLAGS_CPLEX) -I../../src $(INC_MC) $(INC_PROFIL) $(INC_FILIB) $(INC_LAPACK) $(INC_GUROBI) $(INC_CPLEX) $(INC_IPOPT) -o test7.o test7.cpp

test8 : test8.o
	g++ $(PROF) -o test8 test8.o $(LIB_MC) $(LIB_PROFIL) $(LIB_FILIB) $(LIB_HSL) $(LIB_LAPACK) $(LIB_GUROBI) $(LIB_CPLEX) $(LIB_IPOPT)

test8.o: test8.cpp
	g++ -c $(PROF) -g -O3 -Wall -std=c++11 $(FLAGS_FILIB) $(FLAGS_CPLEX) -I../../src $(INC_MC) $(INC_PROFIL) $(INC_FILIB) $(INC_LAPACK) $(INC_GUROBI) $(INC_CPLEX) $(INC_IPOPT) -o test8.o test8.cpp

test9 : test9.o
	g++ $(PROF) -o test9 test9.o $(LIB_MC) $(LIB_PROFIL) $(LIB_FILIB) $(LIB_HSL) $(LIB_LAPACK) $(LIB_GUROBI) $(LIB_IPOPT)

test9.o: test9.cpp
	g++ -c $(PROF) -g -O2 -Wall -std=c++11 $(FLAGS_FILIB) -I../../src $(INC_MC) $(INC_PROFIL) $(INC_FILIB) $(INC_LAPACK) $(INC_GUROBI) $(INC_IPOPT) -o test9.o test9.cpp

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
	rm -f test? test?b *.o *.out *~ *.eps *.lp

