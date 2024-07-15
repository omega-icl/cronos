# THIRD-PARTY LIBRARIES <<-- CHANGE AS APPROPRIATE -->>

PATH_CRONOS = $(shell cd ../.. ; pwd)
include $(PATH_CRONOS)/src/makeoptions.mk

#####

all: test0 test0b test1 test1b test2 test3 test4

test% : test%.o
	$(CPP) $^ $(LIB_CRONOS) $(LIB_DEP) -o $@

%.o : %.cpp
	$(CPP) -c $(FLAG_CPP) $(INC_CRONOS) $(INC_DEP) $< -o $@

#####

graph : gendot gensvg

gendot :
	dot -Tpng -O $(shell ls *.dot)

gensvg :
	@for FILE in $(shell ls *.png); do \
	  echo $$FILE; \
	  convert -resize 70% $$FILE $$FILE; \
	  (display $$FILE &); \
	done

#####

clean:
	rm -f *.o *.out *.dat *.eps *.dot *.svg *.png *~
	rm -f test? test??

