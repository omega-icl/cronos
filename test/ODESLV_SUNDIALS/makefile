# Third-Party Libraries -->> CHANGE ACCORDINGLY <<--
PATH_MC = /home/bchachua/Programs/Devel/MC++/MC++_2.0
LIB_MC = -llapack -lblas
INC_MC = -I$(PATH_MC)/include -I$(PATH_MC)/src/3rdparty/fadbad++ -I$(PATH_MC)/src/3rdparty/cpplapack/include

PATH_SUNDIALS = /opt/sundials-2.5.0
LIB_SUNDIALS = -L$(PATH_SUNDIALS)/lib -lsundials_cvodes -lsundials_nvecserial \
               -llapack -lblas
INC_SUNDIALS = -I$(PATH_SUNDIALS)/include -I$(PATH_SUNDIALS)/include/sundials \
               -I$(PATH_SUNDIALS)/include/cvodes -I$(PATH_SUNDIALS)/include/nvector

#####

all: test1

test1 : test1.o
	g++ -o test1 test1.o $(LIB_MC) $(LIB_SUNDIALS) $(LIB_LAPACK)

test1.o: test1.cpp
	g++ -c -g -Wall -std=c++11 $(FLAGS_FILIB) -I../../src $(INC_MC) $(INC_SUNDIALS) $(INC_LAPACK) -o test1.o test1.cpp

#####

clean :
	rm -f test1 *.o *.out

