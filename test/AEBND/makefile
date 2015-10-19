# Third-Party Libraries -->> CHANGE ACCORDINGLY <<--
PATH_MC = /home/bchachua/Programs/Devel/MC++/MC++_2.0
LIB_MC = -llapack -lblas
#INC_MC = -I$(PATH_MC)/include -I$(PATH_MC)/src/3rdparty/fadbad++ -I$(PATH_MC)/src/3rdparty/cpplapack/include
INC_MC = -I$(PATH_MC)/src/mc -I$(PATH_MC)/src/3rdparty/fadbad++ -I$(PATH_MC)/src/3rdparty/cpplapack/include

LIB_HSL = -lmc13 -lmc21
INC_HSL = 

PATH_PROFIL = /home/bchachua/Programs/ThirdParty/Profil-2.0.8
LIB_PROFIL = -L$(PATH_PROFIL)/lib -lProfilPackages -lProfil -lBias -llr
INC_PROFIL = -I$(PATH_PROFIL)/include

PATH_FILIB = /opt/filib++
LIB_FILIB = -L$(PATH_FILIB)/lib -lprim
INC_FILIB = -I$(PATH_FILIB)/include/ -I$(PATH_FILIB)/include/interval
FLAGS_FILIB = -frounding-math -ffloat-store

#####

all: test1 test2 test3 test4 test5 test6 test7 test8 test9 test10

test1 : test1.o
	g++ -o test1 test1.o $(LIB_MC) $(LIB_PROFIL) $(LIB_FILIB) $(LIB_HSL) $(LIB_LAPACK)

test1.o: test1.cpp
	g++ -c -g -std=c++11 -Wall $(FLAGS_FILIB) -I../../src $(INC_MC) $(INC_PROFIL) $(INC_FILIB) $(INC_LAPACK) -o test1.o test1.cpp

test2 : test2.o
	g++ -o test2 test2.o $(LIB_MC) $(LIB_PROFIL) $(LIB_FILIB) $(LIB_HSL) $(LIB_LAPACK)

test2.o: test2.cpp
	g++ -c -g -std=c++11 -Wall $(FLAGS_FILIB) -I../../src $(INC_MC) $(INC_PROFIL) $(INC_FILIB) $(INC_LAPACK) -o test2.o test2.cpp

test3 : test3.o
	g++ -o test3 test3.o $(LIB_MC) $(LIB_PROFIL) $(LIB_FILIB) $(LIB_HSL) $(LIB_LAPACK)

test3.o: test3.cpp
	g++ -c -g -std=c++11 -Wall $(FLAGS_FILIB) -I../../src $(INC_MC) $(INC_PROFIL) $(INC_FILIB) $(INC_LAPACK) -o test3.o test3.cpp

test4 : test4.o
	g++ -o test4 test4.o $(LIB_MC) $(LIB_PROFIL) $(LIB_FILIB) $(LIB_HSL) $(LIB_LAPACK)

test4.o: test4.cpp
	g++ -c -g -std=c++11 -Wall $(FLAGS_FILIB) -I../../src $(INC_MC) $(INC_PROFIL) $(INC_FILIB) $(INC_LAPACK) -o test4.o test4.cpp

test5 : test5.o
	g++ -o test5 test5.o $(LIB_MC) $(LIB_PROFIL) $(LIB_FILIB) $(LIB_HSL) $(LIB_LAPACK)

test5.o: test5.cpp
	g++ -c -g -std=c++11 -Wall $(FLAGS_FILIB) -I../../src $(INC_MC) $(INC_PROFIL) $(INC_FILIB) $(INC_LAPACK) -o test5.o test5.cpp

test6 : test6.o
	g++ -o test6 test6.o $(LIB_MC) $(LIB_PROFIL) $(LIB_FILIB) $(LIB_HSL) $(LIB_LAPACK)

test6.o: test6.cpp
	g++ -c -g -std=c++11 -O2 -Wall $(FLAGS_FILIB) -I../../src $(INC_MC) $(INC_PROFIL) $(INC_FILIB) $(INC_LAPACK) -o test6.o test6.cpp

test7 : test7.o
	g++ -o test7 test7.o $(LIB_MC) $(LIB_PROFIL) $(LIB_FILIB) $(LIB_HSL) $(LIB_LAPACK)

test7.o: test7.cpp
	g++ -c -g -std=c++11 -O2 -Wall $(FLAGS_FILIB) -I../../src $(INC_MC) $(INC_PROFIL) $(INC_FILIB) $(INC_LAPACK) -o test7.o test7.cpp

test8 : test8.o
	g++ -o test8 test8.o $(LIB_MC) $(LIB_PROFIL) $(LIB_FILIB) $(LIB_HSL) $(LIB_LAPACK)

test8.o: test8.cpp
	g++ -c -g -std=c++11 -Wall $(FLAGS_FILIB) -I../../src $(INC_MC) $(INC_PROFIL) $(INC_FILIB) $(INC_LAPACK) -o test8.o test8.cpp

test9 : test9.o
	g++ -o test9 test9.o $(LIB_MC) $(LIB_PROFIL) $(LIB_FILIB) $(LIB_HSL) $(LIB_LAPACK)

test9.o: test9.cpp
	g++ -c -g -std=c++11 -Wall $(FLAGS_FILIB) -I../../src $(INC_MC) $(INC_PROFIL) $(INC_FILIB) $(INC_LAPACK) -o test9.o test9.cpp

test10 : test10.o
	g++ -o test10 test10.o $(LIB_MC) $(LIB_PROFIL) $(LIB_FILIB) $(LIB_HSL) $(LIB_LAPACK)

test10.o: test10.cpp
	g++ -c -g -std=c++11 -Wall $(FLAGS_FILIB) -I../../src $(INC_MC) $(INC_PROFIL) $(INC_FILIB) $(INC_LAPACK) -o test10.o test10.cpp

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
	rm -f test? test?? *.o *.out *.dot *.svg *~

