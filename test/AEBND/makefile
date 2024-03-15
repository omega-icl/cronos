# THIRD-PARTY LIBRARIES <<-- CHANGE AS APPROPRIATE -->>

PATH_CRONOS = $(shell cd ../.. ; pwd)
LIB_CRONOS = #-L$(PATH_CRONOS)/lib -Wl,-rpath=$(PATH_CRONOS)/lib -lcronos
INC_CRONOS = -I$(PATH_CRONOS)/include

include $(PATH_CRONOS)/src/makeoptions.mk

#####

all: test1 test2 test3 test4 test5 test6

test% : test%.o
	$(CPP) $^ $(LIB_CRONOS) $(LIB_DEP) -o $@

%.o : %.cpp
	$(CPP) -c $(FLAG_CPP) $(INC_CRONOS) $(INC_DEP) $< -o $@

#####

clean:
	rm -f *.o *.out *.dat *.eps *.dot *.svg *~
	rm -f test? test??

