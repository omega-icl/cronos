# Third-Party Libraries -->> CHANGE ACCORDINGLY <<--
PATH_MC = $(HOME)/Programs/Devel/MC++/MC++_2.0
LIB_MC = -llapack -lblas
INC_MC = -I$(PATH_MC)/include -I$(PATH_MC)/src/3rdparty/fadbad++ -I$(PATH_MC)/src/3rdparty/cpplapack/include

LIB_GSL = -lgsl -lgslcblas -lm
INC_GSL = 

PATH_SUNDIALS = /opt/sundials-2.5.0
#PATH_SUNDIALS = /usr/local
LIB_SUNDIALS = -L$(PATH_SUNDIALS)/lib -lsundials_cvodes -lsundials_nvecserial \
               -llapack -lblas
INC_SUNDIALS = -I$(PATH_SUNDIALS)/include -I$(PATH_SUNDIALS)/include/sundials \
               -I$(PATH_SUNDIALS)/include/cvodes -I$(PATH_SUNDIALS)/include/nvector

PATH_PROFIL = $(HOME)/Programs/ThirdParty/Profil-2.0.8
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

all: test1 test2 test3 test4 test5 test5b test6

test% : test%.o
	g++ -o $@ $< $(LIB_ALL)

test%.o: test%.cpp
	g++ -c $(FLAG_ALL) -I../../src $(INC_ALL) -o $@ $<

#####

clean :
	rm -f test? test?b *.o *.out *.dat *.eps *~

