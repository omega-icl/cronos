# Third-Party Libraries -->> CHANGE ACCORDINGLY <<--
PATH_MC = $(HOME)/Programs/Devel/MC++/MC++_2.1
LIB_MC = -llapack -lblas
INC_MC = -I$(PATH_MC)/src/mc -I$(PATH_MC)/src/3rdparty/fadbad++ -I$(PATH_MC)/src/3rdparty/cpplapack-2015.05.11-1/include
FLAG_MC = -Wno-misleading-indentation -Wno-unknown-pragmas -DMC__USE_HSL

PATH_IPOPT   = /opt/CoinIpopt
PATH_COINHSL = /opt/coinhsl
LIB_IPOPT = -L$(PATH_IPOPT)/lib -L$(PATH_COINHSL)/lib -lipopt -lcoinhsl -llapack -lblas -lgfortran -lm -lquadmath -lblas -lgfortran -lm -lquadmath -lm -ldl -lcoinmumps -lblas -lgfortran -lm -lquadmath -lgfortran -lm -lquadmath -lcoinmetis
INC_IPOPT = -I$(PATH_IPOPT)/include
#PATH_IPOPT = /opt/Ipopt-3.10.2
#LIB_IPOPT = -L$(PATH_IPOPT)/lib -lipopt -lcoinhsl -llapack -lblas -ldl
#INC_IPOPT = -I$(PATH_IPOPT)/include

#FLAG_ALL = -g -Wall -std=c++11
FLAG_ALL = -g -Ofast -Wall -std=c++11 $(FLAG_MC)
LIB_ALL  = $(LIB_MC) $(LIB_IPOPT)
INC_ALL  = $(INC_MC) $(INC_IPOPT)

#####

all: test1 test2 test3 test4 test5

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

