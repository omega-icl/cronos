# Third-Party Libraries -->> CHANGE ACCORDINGLY <<--
PATH_MC = $(HOME)/Programs/Devel/MC++/MC++_2.0
LIB_MC = -llapack -lblas
INC_MC = -I$(PATH_MC)/include -I$(PATH_MC)/src/3rdparty/fadbad++ -I$(PATH_MC)/src/3rdparty/cpplapack/include

LIB_HSL = -lmc13 -lmc21 -lmc33 -lgfortran
INC_HSL = 

PATH_PROFIL = /home/bchachua/Programs/ThirdParty/Profil-2.0.8
LIB_PROFIL  = -L$(PATH_PROFIL)/lib -lProfilPackages -lProfil -lBias -llr
INC_PROFIL  = -I$(PATH_PROFIL)/include

PATH_FILIB  = /opt/filib++
LIB_FILIB   = -L$(PATH_FILIB)/lib -lprim
INC_FILIB   = -I$(PATH_FILIB)/include/ -I$(PATH_FILIB)/include/interval
FLAGS_FILIB = -frounding-math -ffloat-store

#PROF = -pg
#OPTI = -Ofast
#OPTI = -O2

FLAG_ALL = -g $(OPTI) -Wall -std=c++11 $(FLAG_FILIB)
LIB_ALL  = $(LIB_MC) $(LIB_PROFIL) $(LIB_FILIB) $(LIB_HSL)
INC_ALL  = $(INC_MC) $(INC_PROFIL) $(INC_FILIB) $(INC_HSL)

#####

all: test1 test2 test3 test4 test5

test% : test%.o
	g++ $(PROF) -o $@ $< $(LIB_ALL)

test%.o: test%.cpp
	g++ $(PROF) -c $(FLAG_ALL) -I../../src $(INC_ALL) -o $@ $<

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

