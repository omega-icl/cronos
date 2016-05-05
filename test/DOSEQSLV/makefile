# Third-Party Libraries -->> CHANGE ACCORDINGLY <<--
PATH_MC = /home/bchachua/Programs/Devel/MC++/MC++_2.0
LIB_MC = -llapack -lblas
INC_MC = -I$(PATH_MC)/include -I$(PATH_MC)/src/3rdparty/fadbad++ -I$(PATH_MC)/src/3rdparty/cpplapack/include

PATH_SUNDIALS = /opt/sundials-2.6.2
LIB_SUNDIALS = -L$(PATH_SUNDIALS)/lib -lsundials_cvodes -lsundials_nvecserial \
               -llapack -lblas
INC_SUNDIALS = -I$(PATH_SUNDIALS)/include -I$(PATH_SUNDIALS)/include/sundials \
               -I$(PATH_SUNDIALS)/include/cvodes -I$(PATH_SUNDIALS)/include/nvector

PATH_IPOPT = /opt/Ipopt-3.10.2
LIB_IPOPT = -L$(PATH_IPOPT)/lib -lipopt -lcoinhsl -llapack -lblas -ldl
INC_IPOPT = -I$(PATH_IPOPT)/include

#FLAG_ALL = -g -Wall -std=c++11
FLAG_ALL = -g -O3 -Wall -std=c++11
LIB_ALL  = $(LIB_MC) $(LIB_IPOPT) $(LIB_SUNDIALS) $(LIB_LAPACK)
INC_ALL  = $(INC_MC) $(INC_IPOPT) $(INC_SUNDIALS) $(INC_LAPACK)

#####

all: test1 test2

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
	rm -f test? test?b *.o *.out *.dat *.eps *.dot *~

