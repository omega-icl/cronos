# Third-Party Libraries -->> CHANGE ACCORDINGLY <<--
PATH_MC = /home/bchachua/Programs/Devel/MC++/MC++_2.0
LIB_MC = -llapack -lblas
INC_MC = -I$(PATH_MC)/include -I$(PATH_MC)/src/3rdparty/fadbad++ -I$(PATH_MC)/src/3rdparty/cpplapack/include

PATH_IPOPT = /opt/Ipopt-3.10.2
LIB_IPOPT = -L$(PATH_IPOPT)/lib -lipopt -lcoinhsl -llapack -lblas -ldl
INC_IPOPT = -I$(PATH_IPOPT)/include

#####

all: test1 test2 test3 test4 test5

test1 : test1.o
	g++ -o test1 test1.o $(LIB_MC) $(LIB_IPOPT)

test1.o: test1.cpp
	g++ -c -g -Wall -std=c++11 -I../../src $(INC_MC) $(INC_IPOPT) -o test1.o test1.cpp

test2 : test2.o
	g++ -o test2 test2.o $(LIB_MC) $(LIB_IPOPT)

test2.o: test2.cpp
	g++ -c -g -Wall -std=c++11 -I../../src $(INC_MC) $(INC_IPOPT) -o test2.o test2.cpp

test3 : test3.o
	g++ -o test3 test3.o $(LIB_MC) $(LIB_IPOPT)

test3.o: test3.cpp
	g++ -c -g -Wall -std=c++11 -I../../src $(INC_MC) $(INC_IPOPT) -o test3.o test3.cpp

test4 : test4.o
	g++ -o test4 test4.o $(LIB_MC) $(LIB_IPOPT)

test4.o: test4.cpp
	g++ -c -g -Wall -std=c++11 -I../../src $(INC_MC) $(INC_IPOPT) -o test4.o test4.cpp

test5 : test5.o
	g++ -o test5 test5.o $(LIB_MC) $(LIB_IPOPT)

test5.o: test5.cpp
	g++ -c -g -Wall -std=c++11 -I../../src $(INC_MC) $(INC_IPOPT) -o test5.o test5.cpp

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
	rm -f test? *.o *.out *.eps *~

