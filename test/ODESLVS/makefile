# Third-Party Libraries -->> CHANGE ACCORDINGLY <<--
PATH_MC = $(HOME)/Programs/Devel/MC++/MC++_2.1
LIB_MC = -llapack -lblas
INC_MC = -I$(PATH_MC)/src/mc -I$(PATH_MC)/src/3rdparty/fadbad++ -I$(PATH_MC)/src/3rdparty/cpplapack-2015.05.11-1/include
FLAG_MC = -Wno-misleading-indentation -Wno-unknown-pragmas

PATH_SUNDIALS = /opt/sundials-2.7.0
LIB_SUNDIALS = -L$(PATH_SUNDIALS)/lib -lsundials_cvodes -lsundials_nvecserial \
               -llapack -lblas
INC_SUNDIALS = -I$(PATH_SUNDIALS)/include -I$(PATH_SUNDIALS)/include/sundials \
               -I$(PATH_SUNDIALS)/include/cvodes -I$(PATH_SUNDIALS)/include/nvector

LIB_GSL = -lgsl -lgslcblas -lm
INC_GSL = 

#FLAG_ALL = -g -Wall -std=c++11 $(FLAG_FILIB)
FLAG_ALL = -g -Ofast -Wall -std=c++11 $(FLAG_FILIB) $(FLAG_MC)
LIB_ALL  = $(LIB_MC) $(LIB_GSL) $(LIB_SUNDIALS) $(LIB_LAPACK)
INC_ALL  = $(INC_MC) $(INC_GSL) $(INC_SUNDIALS) $(INC_LAPACK)

#####

all: test1 test2 test3 test4

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
	rm -f test? test?b *.o *.out *.dat *.eps *.dot *.svg *~

