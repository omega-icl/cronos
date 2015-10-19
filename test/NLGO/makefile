# Third-Party Libraries -->> CHANGE ACCORDINGLY <<--
PATH_MC = /home/bchachua/Programs/Devel/MC++/MC++_2.0
LIB_MC = -llapack -lblas
INC_MC = -I$(PATH_MC)/include -I$(PATH_MC)/src/3rdparty/fadbad++ -I$(PATH_MC)/src/3rdparty/cpplapack/include

PATH_PROFIL = /home/bchachua/Programs/ThirdParty/Profil-2.0.8
LIB_PROFIL  = -L$(PATH_PROFIL)/lib -lProfilPackages -lProfil -lBias -llr
INC_PROFIL  = -I$(PATH_PROFIL)/include

PATH_FILIB  = /opt/filib++
LIB_FILIB   = -L$(PATH_FILIB)/lib -lprim
INC_FILIB   = -I$(PATH_FILIB)/include/ -I$(PATH_FILIB)/include/interval
FLAGS_FILIB = -frounding-math -ffloat-store

PATH_IPOPT = /opt/Ipopt-3.10.2
LIB_IPOPT = -L$(PATH_IPOPT)/lib -lipopt -lcoinhsl -llapack -lblas -ldl
INC_IPOPT = -I$(PATH_IPOPT)/include

PATH_GUROBI = $(GUROBI_HOME)
LIB_GUROBI = -L$(PATH_GUROBI)/lib -lgurobi_c++ -lgurobi56 -pthread
INC_GUROBI = -I$(PATH_GUROBI)/include

#####

all: test1 test2 test3 test4 test5

test1 : test1.o
	g++ -o test1 test1.o $(LIB_PROFIL) $(LIB_FILIB) $(LIB_MC) $(LIB_IPOPT) $(LIB_GUROBI)

test1.o: test1.cpp
	g++ -c -g -Wall -std=c++11 -I../../src $(INC_MC) $(INC_PROFIL) $(INC_FILIB) $(INC_IPOPT) $(INC_GUROBI) -o test1.o test1.cpp

test2 : test2.o
	g++ -o test2 test2.o $(LIB_PROFIL) $(LIB_FILIB) $(LIB_MC) $(LIB_IPOPT) $(LIB_GUROBI)

test2.o: test2.cpp
	g++ -c -g -Wall -std=c++11 -I../../src $(INC_MC) $(INC_PROFIL) $(INC_FILIB) $(INC_IPOPT) $(INC_GUROBI) -o test2.o test2.cpp

test3 : test3.o
	g++ -o test3 test3.o $(LIB_PROFIL) $(LIB_FILIB) $(LIB_MC) $(LIB_IPOPT) $(LIB_GUROBI)

test3.o: test3.cpp
	g++ -c -g -Wall -std=c++11 -I../../src  $(INC_MC) $(INC_PROFIL) $(INC_FILIB) $(INC_IPOPT) $(INC_GUROBI) -o test3.o test3.cpp

test4 : test4.o
	g++ -o test4 test4.o $(LIB_PROFIL) $(LIB_FILIB) $(LIB_MC) $(LIB_IPOPT) $(LIB_GUROBI)

test4.o: test4.cpp
	g++ -c -g -Wall -std=c++11 -I../../src  $(INC_MC) $(INC_PROFIL) $(INC_FILIB) $(INC_IPOPT) $(INC_GUROBI) -o test4.o test4.cpp

test5 : test5.o
	g++ -o test5 test5.o $(LIB_PROFIL) $(LIB_FILIB) $(LIB_MC) $(LIB_IPOPT) $(LIB_GUROBI)

test5.o: test5.cpp
	g++ -c -g -O2 -Wall -std=c++11 -I../../src $(INC_MC) $(INC_PROFIL) $(INC_FILIB) $(INC_IPOPT) $(INC_GUROBI) -o test5.o test5.cpp

test6 : test6.o
	g++ -o test6 test6.o $(LIB_PROFIL) $(LIB_FILIB) $(LIB_MC) $(LIB_IPOPT) $(LIB_GUROBI)

test6.o: test6.cpp
	g++ -c -g -O2 -Wall -std=c++11 -I../../src $(INC_MC) $(INC_PROFIL) $(INC_FILIB) $(INC_IPOPT) $(INC_GUROBI) -o test6.o test6.cpp

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
	rm -f test? *.o *.out *.lp *.log *~

