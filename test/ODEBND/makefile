# Third-Party Libraries -->> CHANGE ACCORDINGLY <<--
PATH_MC = /home/bchachua/Programs/Devel/MC++/MC++_2.0
LIB_MC = -llapack -lblas
INC_MC = -I$(PATH_MC)/include -I$(PATH_MC)/src/3rdparty/fadbad++ -I$(PATH_MC)/src/3rdparty/cpplapack/include

LIB_HSL = -lmc13 -lmc21
INC_HSL = 

LIB_GSL = -lgsl -lgslcblas -lm
INC_GSL = 

PATH_SUNDIALS = /opt/sundials-2.5.0
#PATH_SUNDIALS = /usr/local
LIB_SUNDIALS = -L$(PATH_SUNDIALS)/lib -lsundials_cvodes -lsundials_nvecserial \
               -llapack -lblas
INC_SUNDIALS = -I$(PATH_SUNDIALS)/include -I$(PATH_SUNDIALS)/include/sundials \
               -I$(PATH_SUNDIALS)/include/cvodes -I$(PATH_SUNDIALS)/include/nvector

PATH_PROFIL = /home/bchachua/Programs/ThirdParty/Profil-2.0.8
LIB_PROFIL = -L$(PATH_PROFIL)/lib -lProfilPackages -lProfil -lBias -llr
INC_PROFIL = -I$(PATH_PROFIL)/include

PATH_FILIB = /opt/filib++
LIB_FILIB  = -L$(PATH_FILIB)/lib -lprim
INC_FILIB  = -I$(PATH_FILIB)/include/ -I$(PATH_FILIB)/include/interval
FLAG_FILIB = -frounding-math -ffloat-store

FLAG_ALL = -g -Wall -std=c++11 $(FLAG_FILIB)
LIB_ALL  = $(LIB_MC) $(LIB_PROFIL) $(LIB_FILIB) $(LIB_GSL) $(LIB_SUNDIALS) $(LIB_LAPACK)
INC_ALL  = $(INC_MC) $(INC_PROFIL) $(INC_FILIB) $(INC_GSL) $(INC_SUNDIALS) $(INC_LAPACK)

#####

all: test0 test1 test2 test3 test4 test5

test0 : test0.o
	g++ -o test0 test0.o $(LIB_ALL)
 
test0.o: test0.cpp
	g++ -c $(FLAG_ALL) -I../../src $(INC_ALL) -o test0.o test0.cpp

test1 : test1.o
	g++ -o test1 test1.o $(LIB_ALL)
 
test1.o: test1.cpp
	g++ -c $(FLAG_ALL) -I../../src $(INC_ALL) -o test1.o test1.cpp

test2 : test2.o
	g++ -o test2 test2.o $(LIB_ALL)
 
test2.o: test2.cpp
	g++ -c $(FLAG_ALL) -I../../src $(INC_ALL) -o test2.o test2.cpp

test3 : test3.o
	g++ -o test3 test3.o $(LIB_ALL)
 
test3.o: test3.cpp
	g++ -c $(FLAG_ALL) -I../../src $(INC_ALL) -o test3.o test3.cpp

test4 : test4.o
	g++ -o test4 test4.o $(LIB_ALL)
 
test4.o: test4.cpp
	g++ -c $(FLAG_ALL) -I../../src $(INC_ALL) -o test4.o test4.cpp

test5 : test5.o
	g++ -o test5 test5.o $(LIB_ALL)
 
test5.o: test5.cpp
	g++ -c $(FLAG_ALL) -I../../src $(INC_ALL) -o test5.o test5.cpp

#####

clean :
	rm -f test0 test1 test2 test3 test4 test5 *.o *.out *.dat *.eps *~

