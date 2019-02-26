# This makefile compiles a chared library of CRONOS and it creates symbolic links
# to the library in $(libpath) and to the header files in $(incpath)

include $(srcpath)/makeoptions.mk

#####

incobjs = base_ae.hpp base_opt.hpp base_nlp.hpp base_de.hpp base_do.hpp base_rk.hpp \
          base_expand.hpp base_sundials.hpp aebnd.hpp odeslv_base.hpp odeslv_sundials.hpp \
          odeslvs_base.hpp odeslvs_sundials.hpp odebnd_base.hpp odebnd_sundials.hpp \
          odebnd_val.hpp odebnd_expand.hpp iodebnd_base.hpp iodebnd_sundials.hpp \
          odebnds_base.hpp odebnds_sundials.hpp nlpslv_ipopt.hpp doseqslv_ipopt.hpp \
          sbb.hpp sbp.hpp lprelax_base.hpp csearch_base.hpp nlgo.hpp nlcp.hpp doseqgo.hpp

libobjs = odeslvs_sundials.o odeslv_sundials.o
#libobjs = base_de.o base_sundials.o odeslv_base.o odeslv_sundials.o odeslvs_base.o \
#          odeslvs_sundials.o

libname = libcronos.so

#####

install: dispBuild sobol_lib cronos_lib dispInstall
	@if test ! -e $(libpath)/$(libname); then \
		echo creating symolic link to shared library $(libname); \
		cd $(libpath) ; ln -s $(srcpath)/$(libname) $(libname); \
	fi
	@for INC in $(incobjs); do \
		if test ! -e $(incpath)/$$INC; then \
			echo creating symbolic link to header file $$INC; \
			cd $(incpath); ln -s $(srcpath)/$$INC $$INC; \
		fi; \
	done
	@echo

sobol_lib:
	-(cd $(thirdppath); make -f sobol.mk lib)

cronos_lib: $(libobjs)
	$(CPP) -shared -o $(libname) $(libobjs) -L$(thirdppath) -lsobol

%.o : %.cpp
	$(CPP) -c $(FLAG_CPP) $(INC_DEP) $< -o $@

dispBuild:
	@echo
	@(echo '***Compiling CRONOS library (ver.' $(version)')***')
	@echo

dispInstall:
	@echo
	@(echo '***Installing CRONOS library (ver.' $(version)')***')
	@echo

#####

clean: dispClean
	rm -fi $(libobjs) $(libname)

dispClean:
	@echo
	@(echo '***Cleaning CRONOS directory (ver.' $(version)')***')
	@echo

#####

cleandist: dispCleanInstall
	cd $(thirdppath); make -f sobol.mk clean
	rm -fi $(libobjs) $(libname)
	cd $(incpath) ; rm -i $(incobjs)
	cd $(libpath) ; rm -i $(libname)

dispCleanInstall:
	@echo
	@(echo '***Uninstalling CRONOS library (ver.' $(version)')***')
	@echo
