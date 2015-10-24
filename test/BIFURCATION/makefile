# Third-Party Libraries -->> CHANGE ACCORDINGLY <<--
PATH_MC = /home/bchachua/Programs/Devel/MC++/MC++_2.0
LIB_MC = -llapack -lblas
INC_MC = -I$(PATH_MC)/include -I$(PATH_MC)/src/3rdparty/fadbad++ -I$(PATH_MC)/src/3rdparty/cpplapack/include

LIB_HSL = -lmc13
INC_HSL = 

PATH_PROFIL = /home/bchachua/Programs/ThirdParty/Profil-2.0.8
LIB_PROFIL  = -L$(PATH_PROFIL)/lib -lProfilPackages -lProfil -lBias -llr
INC_PROFIL  = -I$(PATH_PROFIL)/include

PATH_FILIB  = /opt/filib++
LIB_FILIB   = -L$(PATH_FILIB)/lib -lprim
INC_FILIB   = -I$(PATH_FILIB)/include/ -I$(PATH_FILIB)/include/interval
FLAGS_FILIB = -frounding-math -ffloat-store

PATH_GUROBI = $(GUROBI_HOME)
LIB_GUROBI  = -L$(PATH_GUROBI)/lib -lgurobi_c++ -lgurobi56 -pthread
INC_GUROBI  = -I$(PATH_GUROBI)/include

LIB_SDPA   = -lsdpa -ldmumps_seq
INC_SDPA   = 

#####

all: test1 test2 test3

test1 : test1.o
	g++ -o test1 test1.o $(LIB_MC) $(LIB_PROFIL) $(LIB_FILIB) $(LIB_HSL) $(LIB_LAPACK) $(LIB_GUROBI) $(LIB_SDPA)

test1.o: test1.cpp
	g++ -c -g -O2 -Wall -std=c++11 $(FLAGS_FILIB) -I../../src $(INC_MC) $(INC_PROFIL) $(INC_FILIB) $(INC_LAPACK) $(INC_GUROBI) $(INC_SDPA) -o test1.o test1.cpp

test2 : test2.o
	g++ -o test2 test2.o $(LIB_MC) $(LIB_PROFIL) $(LIB_FILIB) $(LIB_HSL) $(LIB_LAPACK) $(LIB_GUROBI) $(LIB_SDPA)

test2.o: test2.cpp
	g++ -c -g -O2 -Wall -std=c++11 $(FLAGS_FILIB) -I../../src $(INC_MC) $(INC_PROFIL) $(INC_FILIB) $(INC_LAPACK) $(INC_GUROBI) $(INC_SDPA) -o test2.o test2.cpp

test2b : test2b.o
	g++ -o test2b test2b.o $(LIB_MC) $(LIB_PROFIL) $(LIB_FILIB) $(LIB_HSL) $(LIB_LAPACK) $(LIB_GUROBI) $(LIB_SDPA)

test2b.o: test2b.cpp
	g++ -c -g -O2 -Wall -std=c++11 $(FLAGS_FILIB) -I../../src $(INC_MC) $(INC_PROFIL) $(INC_FILIB) $(INC_LAPACK) $(INC_GUROBI) $(INC_SDPA) -o test2b.o test2b.cpp

test2c : test2c.o
	g++ -o test2c test2c.o $(LIB_MC) $(LIB_PROFIL) $(LIB_FILIB) $(LIB_HSL) $(LIB_LAPACK) $(LIB_GUROBI) $(LIB_SDPA)

test2c.o: test2c.cpp
	g++ -c -g -O2 -Wall -std=c++11 $(FLAGS_FILIB) -I../../src $(INC_MC) $(INC_PROFIL) $(INC_FILIB) $(INC_LAPACK) $(INC_GUROBI) $(INC_SDPA) -o test2c.o test2c.cpp

test2d : test2d.o
	g++ -o test2d test2d.o $(LIB_MC) $(LIB_PROFIL) $(LIB_FILIB) $(LIB_HSL) $(LIB_LAPACK) $(LIB_GUROBI) $(LIB_SDPA)

test2d.o: test2d.cpp
	g++ -c -g -O2 -Wall -std=c++11 $(FLAGS_FILIB) -I../../src $(INC_MC) $(INC_PROFIL) $(INC_FILIB) $(INC_LAPACK) $(INC_GUROBI) $(INC_SDPA) -o test2d.o test2d.cpp

test3 : test3.o
	g++ -o test3 test3.o $(LIB_MC) $(LIB_PROFIL) $(LIB_FILIB) $(LIB_HSL) $(LIB_LAPACK) $(LIB_GUROBI) $(LIB_SDPA)

test3.o: test3.cpp
	g++ -c -g -O2 -Wall -std=c++11 $(FLAGS_FILIB) -I../../src $(INC_MC) $(INC_PROFIL) $(INC_FILIB) $(INC_LAPACK) $(INC_GUROBI) $(INC_SDPA) -o test3.o test3.cpp

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
	rm -f test? *.o *.out

