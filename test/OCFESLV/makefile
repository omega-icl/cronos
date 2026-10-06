PATH_CRONOS = $(shell cd ../.. ; pwd)

# Compilation options
include $(PATH_CRONOS)/src/makeoptions.mk

#####

all: OCFE_revision OCFE_options OCFE_nonuniform OCFE_refinit OCFE_fadbad OCFE_deepcopy OCFE_ffocfe \
     OCFE_parallel OCFE_marching OCFE_scalar OCFE_autoelim OCFE_partialfold OCFE_nlpartial \
     OCFE_strongcorner OCFE_hypclosure OCFE_homotopy OCFE_homotopy_reset OCFE_homotopy_copy \
     OCFE_bisect OCFE_scalar_Tlo OCFE_scalar_point OCFE_scalar_opeval OCFE_scalar_mixed \
     OCFE_scalar_fctopeval OCFE_scalar_nest OCFE_scalar_nest2 OCFE_scalar_nest3 OCFE_scalar_nest4 \
     OCFE_march_capture OCFE_march_capture2 OCFE_march_capture3 OCFE_march_nested OCFE_march_bolza \
     OCFE_mono_capture OCFE_mono_capture2 OCFE_mono_eval_direct OCFE_capture_nested \
     OCFE_capture_distributed OCFE_capture_distributed2 OCFE_capture_mixed OCFE_capture_expr \
     OCFE_capture_marching OCFE_hybrid_adequacy OCFE_dagcopy_dump OCFE_exceptions_api OCFE_resetup \
     OCFE_inpcolloc OCFE_claim OCFE_closure OCFE_receiver OCFE_ffmodel OCFE_setmodel OCFE_zindex \
     OCFE_faces OCFE_pencil OCFE_fixinput OCFE_symdiff1 OCFE_symdiff2 OCFE_evalsum OCFE_ffocfe2 \
     OCFE_TRANSoc OCFE_TRANSpde OCFE_TRANSpsa OCFE_TRANSdae OCFE_TRANSode OCFE_TRANSgrad \
     OCFE_TRANSeval OCFE_MIXEDred OCFE_AUTOrole OCFE_causal OCFE_marchflux \
     OCFE_AE0 OCFE_AE1 OCFE_AE2 \
     OCFE_ODE0 OCFE_ODE1 OCFE_ODE2 \
     OCFE_DAE0 OCFE_DAE1 OCFE_DAE2 OCFE_DAE3 OCFE_DAE4 OCFE_DAE5 OCFE_DAE5sd \
     OCFE_ALGFIELD1 OCFE_ALGFIELD2 \
     OCFE_PDE0 OCFE_PDE1 OCFE_PDE2 OCFE_PDE3 OCFE_PDE4 OCFE_PDE5 OCFE_PDE5b OCFE_PDE6 OCFE_PDE6b \
     OCFE_PDE7 OCFE_PDE8 OCFE_PDE9 OCFE_PDE10 OCFE_PDE11 OCFE_PDE12 OCFE_PDE12m OCFE_PDE13 OCFE_PDE14 \
     OCFE_PDE15 OCFE_PDE16 OCFE_PDE17 OCFE_PDE18 OCFE_PDE19 OCFE_PDE20 OCFE_PDE20b OCFE_PDE20c \
     OCFE_PDE20d OCFE_PDE20e OCFE_PDE20f OCFE_PDE20g OCFE_PDE21 OCFE_PDE22 OCFE_PDE23 OCFE_PDE23b \
     OCFE_PDE24 OCFE_PDE25 OCFE_PDE26 OCFE_PDE27 OCFE_PDE28 OCFE_PDE29 OCFE_PDE30 OCFE_PDE30b \
     OCFE_PDE31 OCFE_PDE31b OCFE_PDE32 OCFE_PDE32b OCFE_PDE32c OCFE_PDE33 OCFE_PDE34 OCFE_PDE35 \
     OCFE_PDE36 OCFE_PDE36b OCFE_PDE36c OCFE_PDE37 OCFE_PDE38 OCFE_PDE39 OCFE_PDE39b OCFE_PDE40 \
     OCFE_PDE40b \
     OCFE_CSTR1 OCFE_CSTR2a OCFE_CSTR2b OCFE_CSTR3 OCFE_CSTR4 OCFE_CSTR4b \
     OCFE_PSA1 OCFE_PSA2 OCFE_PSA3 OCFE_PSA4 OCFE_PSA5 OCFE_PSA6 OCFE_PSA7 OCFE_PSA8 OCFE_PSA9 \
     OCFE_PSA10 \
     OCFE_MBC1 OCFE_MBC2a OCFE_MBC2b OCFE_MBC3 OCFE_MBC4 OCFE_MBC5 OCFE_MBC5b OCFE_MBC5c \
     OCFE_MBC5fig OCFE_MBC5fig_ebal OCFE_MBC6 OCFE_MBC6fig OCFE_MBC6fig_ebal \
     OCFE_MMPDE1 OCFE_MESHONLY1

OCFE_revision : OCFE_revision.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_options : OCFE_options.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_nonuniform : OCFE_nonuniform.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_refinit : OCFE_refinit.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_fadbad : OCFE_fadbad.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_deepcopy : OCFE_deepcopy.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_ffocfe : OCFE_ffocfe.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_parallel : OCFE_parallel.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_marching : OCFE_marching.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_scalar : OCFE_scalar.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_autoelim : OCFE_autoelim.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_partialfold : OCFE_partialfold.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_nlpartial : OCFE_nlpartial.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_strongcorner : OCFE_strongcorner.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_hypclosure : OCFE_hypclosure.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_homotopy : OCFE_homotopy.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_homotopy_reset : OCFE_homotopy_reset.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_homotopy_copy : OCFE_homotopy_copy.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_bisect : OCFE_bisect.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)


OCFE_scalar_Tlo : OCFE_scalar_Tlo.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_scalar_point : OCFE_scalar_point.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_scalar_opeval : OCFE_scalar_opeval.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_scalar_mixed : OCFE_scalar_mixed.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_scalar_fctopeval : OCFE_scalar_fctopeval.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_scalar_nest : OCFE_scalar_nest.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_scalar_nest2 : OCFE_scalar_nest2.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_scalar_nest3 : OCFE_scalar_nest3.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_scalar_nest4 : OCFE_scalar_nest4.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)


OCFE_march_capture : OCFE_march_capture.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_march_capture2 : OCFE_march_capture2.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_march_capture3 : OCFE_march_capture3.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_march_nested : OCFE_march_nested.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_march_bolza : OCFE_march_bolza.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_mono_capture : OCFE_mono_capture.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_mono_capture2 : OCFE_mono_capture2.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_mono_eval_direct : OCFE_mono_eval_direct.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_capture_nested : OCFE_capture_nested.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_capture_distributed : OCFE_capture_distributed.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_capture_distributed2 : OCFE_capture_distributed2.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_capture_mixed : OCFE_capture_mixed.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_capture_expr : OCFE_capture_expr.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_capture_marching : OCFE_capture_marching.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_hybrid_adequacy : OCFE_hybrid_adequacy.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)


OCFE_dagcopy_dump : OCFE_dagcopy_dump.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_exceptions_api : OCFE_exceptions_api.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_resetup : OCFE_resetup.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_inpcolloc : OCFE_inpcolloc.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_claim : OCFE_claim.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_closure : OCFE_closure.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_receiver : OCFE_receiver.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_ffmodel : OCFE_ffmodel.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_setmodel : OCFE_setmodel.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_zindex : OCFE_zindex.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_faces : OCFE_faces.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_pencil : OCFE_pencil.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_fixinput : OCFE_fixinput.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_symdiff1 : OCFE_symdiff1.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_symdiff2 : OCFE_symdiff2.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_evalsum : OCFE_evalsum.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_ffocfe2 : OCFE_ffocfe2.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_TRANSoc : OCFE_TRANSoc.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_TRANSpde : OCFE_TRANSpde.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_TRANSpsa : OCFE_TRANSpsa.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_TRANSdae : OCFE_TRANSdae.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_TRANSode : OCFE_TRANSode.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_TRANSgrad : OCFE_TRANSgrad.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_TRANSeval : OCFE_TRANSeval.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_MIXEDred : OCFE_MIXEDred.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_AUTOrole : OCFE_AUTOrole.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_causal : OCFE_causal.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_marchflux : OCFE_marchflux.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)


OCFE_AE0 : OCFE_AE0.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_AE1 : OCFE_AE1.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_AE2 : OCFE_AE2.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)


OCFE_ODE0 : OCFE_ODE0.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_ODE1 : OCFE_ODE1.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_ODE2 : OCFE_ODE2.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)


OCFE_DAE0 : OCFE_DAE0.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_DAE1 : OCFE_DAE1.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_DAE2 : OCFE_DAE2.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_DAE3 : OCFE_DAE3.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_DAE4 : OCFE_DAE4.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_DAE5 : OCFE_DAE5.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_DAE5sd : OCFE_DAE5sd.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)


OCFE_ALGFIELD1 : OCFE_ALGFIELD1.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_ALGFIELD2 : OCFE_ALGFIELD2.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)


OCFE_PDE0 : OCFE_PDE0.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE1 : OCFE_PDE1.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE2 : OCFE_PDE2.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE3 : OCFE_PDE3.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE4 : OCFE_PDE4.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE5 : OCFE_PDE5.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE5b : OCFE_PDE5b.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE6 : OCFE_PDE6.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE6b : OCFE_PDE6b.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE7 : OCFE_PDE7.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE8 : OCFE_PDE8.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE9 : OCFE_PDE9.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE10 : OCFE_PDE10.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE11 : OCFE_PDE11.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE12 : OCFE_PDE12.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE12m : OCFE_PDE12m.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE13 : OCFE_PDE13.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE14 : OCFE_PDE14.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE15 : OCFE_PDE15.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE16 : OCFE_PDE16.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE17 : OCFE_PDE17.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE18 : OCFE_PDE18.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE19 : OCFE_PDE19.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE20 : OCFE_PDE20.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE20b : OCFE_PDE20b.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE20c : OCFE_PDE20c.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE20d : OCFE_PDE20d.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE20e : OCFE_PDE20e.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE20f : OCFE_PDE20f.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE20g : OCFE_PDE20g.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE21 : OCFE_PDE21.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE22 : OCFE_PDE22.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE23 : OCFE_PDE23.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE23b : OCFE_PDE23b.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE24 : OCFE_PDE24.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE25 : OCFE_PDE25.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE26 : OCFE_PDE26.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE27 : OCFE_PDE27.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE28 : OCFE_PDE28.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE29 : OCFE_PDE29.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE30 : OCFE_PDE30.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE30b : OCFE_PDE30b.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE31 : OCFE_PDE31.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE31b : OCFE_PDE31b.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE32 : OCFE_PDE32.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE32b : OCFE_PDE32b.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE32c : OCFE_PDE32c.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE33 : OCFE_PDE33.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE34 : OCFE_PDE34.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE35 : OCFE_PDE35.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE36 : OCFE_PDE36.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE36b : OCFE_PDE36b.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE36c : OCFE_PDE36c.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE37 : OCFE_PDE37.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE38 : OCFE_PDE38.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE39 : OCFE_PDE39.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE39b : OCFE_PDE39b.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE40 : OCFE_PDE40.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PDE40b : OCFE_PDE40b.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)


OCFE_CSTR1 : OCFE_CSTR1.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_CSTR2a : OCFE_CSTR2a.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_CSTR2b : OCFE_CSTR2b.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_CSTR3 : OCFE_CSTR3.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_CSTR4 : OCFE_CSTR4.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_CSTR4b : OCFE_CSTR4b.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)


OCFE_PSA1 : OCFE_PSA1.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PSA2 : OCFE_PSA2.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PSA3 : OCFE_PSA3.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PSA4 : OCFE_PSA4.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PSA5 : OCFE_PSA5.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PSA6 : OCFE_PSA6.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PSA7 : OCFE_PSA7.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PSA8 : OCFE_PSA8.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PSA9 : OCFE_PSA9.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_PSA10 : OCFE_PSA10.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)


OCFE_MBC1 : OCFE_MBC1.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_MBC2a : OCFE_MBC2a.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_MBC2b : OCFE_MBC2b.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_MBC3 : OCFE_MBC3.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_MBC4 : OCFE_MBC4.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_MBC5 : OCFE_MBC5.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_MBC5b : OCFE_MBC5b.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_MBC5c : OCFE_MBC5c.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_MBC5fig : OCFE_MBC5fig.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_MBC5fig_ebal : OCFE_MBC5fig_ebal.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_MBC6 : OCFE_MBC6.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_MBC6fig : OCFE_MBC6fig.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_MBC6fig_ebal : OCFE_MBC6fig_ebal.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)


OCFE_MMPDE1 : OCFE_MMPDE1.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)

OCFE_MESHONLY1 : OCFE_MESHONLY1.o
	$(LINK) $(FLAG_LINK) -o $@ $< $(LIB_CRONOS) $(LIB_DEP)


%.o: %.cpp
#	$(CPP) -c $(FLAG_CPP) $(FLAG_CRONOS) $(INC_CRONOS) $(INC_DEP) $< -o $@
	$(CPP) -c $(FLAG_CPP) -DOCFE_OCFESLV_HEADER='"ocfeslv.hpp"' $(FLAG_CRONOS) $(INC_CRONOS) $(INC_DEP) $< -o $@

#####

clean :
	rm -f  OCFE_revision OCFE_options OCFE_nonuniform OCFE_refinit OCFE_fadbad OCFE_deepcopy \
	OCFE_ffocfe OCFE_parallel OCFE_marching OCFE_scalar OCFE_autoelim OCFE_partialfold OCFE_nlpartial \
	OCFE_strongcorner OCFE_hypclosure OCFE_homotopy OCFE_homotopy_reset OCFE_homotopy_copy \
	OCFE_bisect OCFE_scalar_Tlo OCFE_scalar_point OCFE_scalar_opeval OCFE_scalar_mixed \
	OCFE_scalar_fctopeval OCFE_scalar_nest OCFE_scalar_nest2 OCFE_scalar_nest3 OCFE_scalar_nest4 \
	OCFE_march_capture OCFE_march_capture2 OCFE_march_capture3 OCFE_march_nested OCFE_march_bolza \
	OCFE_mono_capture OCFE_mono_capture2 OCFE_mono_eval_direct OCFE_capture_nested \
	OCFE_capture_distributed OCFE_capture_distributed2 OCFE_capture_mixed OCFE_capture_expr \
	OCFE_capture_marching OCFE_hybrid_adequacy OCFE_dagcopy_dump OCFE_exceptions_api OCFE_resetup \
	OCFE_inpcolloc OCFE_claim OCFE_closure OCFE_receiver OCFE_ffmodel OCFE_setmodel OCFE_zindex \
	OCFE_faces OCFE_pencil OCFE_fixinput OCFE_symdiff1 OCFE_symdiff2 OCFE_evalsum OCFE_ffocfe2 \
	OCFE_TRANSoc OCFE_TRANSpde OCFE_TRANSpsa OCFE_TRANSdae OCFE_TRANSode OCFE_TRANSgrad \
	OCFE_TRANSeval OCFE_MIXEDred OCFE_AUTOrole OCFE_causal OCFE_marchflux \
	OCFE_AE0 OCFE_AE1 OCFE_AE2 \
	OCFE_ODE0 OCFE_ODE1 OCFE_ODE2 \
	OCFE_DAE0 OCFE_DAE1 OCFE_DAE2 OCFE_DAE3 OCFE_DAE4 OCFE_DAE5 OCFE_DAE5sd \
	OCFE_ALGFIELD1 OCFE_ALGFIELD2 \
	OCFE_PDE0 OCFE_PDE1 OCFE_PDE2 OCFE_PDE3 OCFE_PDE4 OCFE_PDE5 OCFE_PDE5b OCFE_PDE6 OCFE_PDE6b \
	OCFE_PDE7 OCFE_PDE8 OCFE_PDE9 OCFE_PDE10 OCFE_PDE11 OCFE_PDE12 OCFE_PDE12m OCFE_PDE13 OCFE_PDE14 \
	OCFE_PDE15 OCFE_PDE16 OCFE_PDE17 OCFE_PDE18 OCFE_PDE19 OCFE_PDE20 OCFE_PDE20b OCFE_PDE20c \
	OCFE_PDE20d OCFE_PDE20e OCFE_PDE20f OCFE_PDE20g OCFE_PDE21 OCFE_PDE22 OCFE_PDE23 OCFE_PDE23b \
	OCFE_PDE24 OCFE_PDE25 OCFE_PDE26 OCFE_PDE27 OCFE_PDE28 OCFE_PDE29 OCFE_PDE30 OCFE_PDE30b \
	OCFE_PDE31 OCFE_PDE31b OCFE_PDE32 OCFE_PDE32b OCFE_PDE32c OCFE_PDE33 OCFE_PDE34 OCFE_PDE35 \
	OCFE_PDE36 OCFE_PDE36b OCFE_PDE36c OCFE_PDE37 OCFE_PDE38 OCFE_PDE39 OCFE_PDE39b OCFE_PDE40 \
	OCFE_PDE40b \
	OCFE_CSTR1 OCFE_CSTR2a OCFE_CSTR2b OCFE_CSTR3 OCFE_CSTR4 OCFE_CSTR4b \
	OCFE_PSA1 OCFE_PSA2 OCFE_PSA3 OCFE_PSA4 OCFE_PSA5 OCFE_PSA6 OCFE_PSA7 OCFE_PSA8 OCFE_PSA9 \
	OCFE_PSA10 \
	OCFE_MBC1 OCFE_MBC2a OCFE_MBC2b OCFE_MBC3 OCFE_MBC4 OCFE_MBC5 OCFE_MBC5b OCFE_MBC5c \
	OCFE_MBC5fig OCFE_MBC5fig_ebal OCFE_MBC6 OCFE_MBC6fig OCFE_MBC6fig_ebal \
	OCFE_MMPDE1 OCFE_MESHONLY1
	rm -f *.o *.out *.dat *.eps *.dot *.png *.lp *.gp *.log *~

