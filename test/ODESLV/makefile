PATH_CRONOS = $(shell cd ../.. ; pwd)

# Compilation options
include $(PATH_CRONOS)/src/makeoptions.mk

#####

all: ODESLV_test1_model_only ODESLV_test2_model_only ODESLV_test1 ODESLV_test2 ODESLV_test3 ODESLV_test4 \
     ODESLV_test5 ODESLV_test6 ODESLV_test7 ODESLV_test8 ODESLV_test9 ODESLV_test10 ODESLV_test11 ODESLV_test12 \
     ODESLV_test14 ODESLV_showsens ODESLV_nonuniform ODESLV_evalsum ODESLV_LVthreads ODESLV_NIFTE \
     ODESLV_massmatrix ODESLV_TRANSode ODESLV_TRANSmulti ODESLV_TRANSeval ODESLV_TRANSjumps ODESLV_EVALadjoint \
     ODESLV_INPUTSgate ODESLV_OPSgate ODESLV_SIDESgate ODESLV_FFODElt


ODESLV_test1_model_only : ODESLV_test1_model_only.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

ODESLV_test2_model_only : ODESLV_test2_model_only.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

ODESLV_test1 : ODESLV_test1.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

ODESLV_test2 : ODESLV_test2.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

ODESLV_test3 : ODESLV_test3.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

ODESLV_test4 : ODESLV_test4.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

ODESLV_test5 : ODESLV_test5.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

ODESLV_test6 : ODESLV_test6.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

ODESLV_test7 : ODESLV_test7.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

ODESLV_test8 : ODESLV_test8.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

ODESLV_test9 : ODESLV_test9.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

ODESLV_test10 : ODESLV_test10.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

ODESLV_test11 : ODESLV_test11.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

ODESLV_test12 : ODESLV_test12.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

ODESLV_test14 : ODESLV_test14.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

ODESLV_showsens : ODESLV_showsens.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

ODESLV_massmatrix : ODESLV_massmatrix.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

ODESLV_nonuniform : ODESLV_nonuniform.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

ODESLV_evalsum : ODESLV_evalsum.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

ODESLV_LVthreads : ODESLV_LVthreads.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

ODESLV_NIFTE : ODESLV_NIFTE.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

ODESLV_TRANSode : ODESLV_TRANSode.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

ODESLV_TRANSmulti : ODESLV_TRANSmulti.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

ODESLV_TRANSeval : ODESLV_TRANSeval.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

ODESLV_TRANSjumps : ODESLV_TRANSjumps.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

ODESLV_EVALadjoint : ODESLV_EVALadjoint.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

ODESLV_INPUTSgate : ODESLV_INPUTSgate.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

ODESLV_OPSgate : ODESLV_OPSgate.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

ODESLV_SIDESgate : ODESLV_SIDESgate.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

ODESLV_FFODElt : ODESLV_FFODElt.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)


%.o : %.cpp
	$(CPP) -c $(FLAG_CPP) $(FLAG_CRONOS) $(INC_CRONOS) $(INC_DEP) $< -o $@

#####

clean :
	rm -f ODESLV_test1_model_only ODESLV_test2_model_only ODESLV_test1 ODESLV_test2 ODESLV_test3 ODESLV_test4 \
	      ODESLV_test5 ODESLV_test6 ODESLV_test7 ODESLV_test8 ODESLV_test9 ODESLV_test10 ODESLV_test11 ODESLV_test12 \
	      ODESLV_test14 ODESLV_showsens ODESLV_nonuniform ODESLV_evalsum ODESLV_LVthreads ODESLV_NIFTE \
	      ODESLV_massmatrix ODESLV_TRANSode ODESLV_TRANSmulti ODESLV_TRANSeval ODESLV_TRANSjumps ODESLV_EVALadjoint \
	      ODESLV_INPUTSgate ODESLV_OPSgate ODESLV_SIDESgate ODESLV_FFODElt
	rm -f *.o *.out *.dat *.eps *.dot *.png *.lp *.gp *.log *~

