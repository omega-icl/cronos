PATH_CRONOS = $(shell cd ../.. ; pwd)

# Compilation options
include $(PATH_CRONOS)/src/makeoptions.mk

#####

all: DAESLV_test1 DAESLV_test2 DAESLV_switched DAESLV_jacobian DAESLV_inputs DAESLV_outputs DAESLV_restart DAESLV_hidden DAESLV_icseed DAESLV_trans DAESLV_icsens DAESLV_fsens

DAESLV_test1 : DAESLV_test1.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

DAESLV_test2 : DAESLV_test2.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

DAESLV_switched : DAESLV_switched.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

DAESLV_jacobian : DAESLV_jacobian.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

DAESLV_inputs : DAESLV_inputs.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

DAESLV_outputs : DAESLV_outputs.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

DAESLV_restart : DAESLV_restart.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

DAESLV_hidden : DAESLV_hidden.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

DAESLV_icseed : DAESLV_icseed.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

DAESLV_trans : DAESLV_trans.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

DAESLV_icsens : DAESLV_icsens.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

DAESLV_fsens : DAESLV_fsens.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

%.o : %.cpp
	$(CPP) -c $(FLAG_CPP) $(FLAG_CRONOS) $(INC_CRONOS) $(INC_DEP) $< -o $@

#####

clean :
	rm -f DAESLV_test1 DAESLV_test2 DAESLV_switched DAESLV_jacobian DAESLV_inputs DAESLV_outputs DAESLV_restart DAESLV_hidden DAESLV_icseed DAESLV_trans DAESLV_icsens DAESLV_fsens
	rm -f *.o *.out *.dat *.eps *.dot *.png *.lp *.gp *.log *~
