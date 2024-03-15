# This makefile compiles a chared library of CRONOS and it creates symbolic links
# to the library in $(libpath) and to the header files in $(incpath)

include $(srcpath)/makeoptions.mk

#####

incobjs = base_ae.hpp base_de.hpp base_rk.hpp base_expand.hpp base_sundials.hpp base_cvodes.hpp \
          aebnd.hpp odeslv_base.hpp odeslv_cvodes.hpp odeslvs_base.hpp odeslvs_cvodes.hpp \
          odebnd_base.hpp odebnd_sundials.hpp odebnd_val.hpp odebnd_expand.hpp \
          iodebnd_base.hpp iodebnd_sundials.hpp odebnds_base.hpp odebnds_sundials.hpp

libobjs =
#libobjs = odeslvs_sundials.o odeslv_sundials.o

libname =
#libname = libcronos.so

#####

install: dispBuild cronos_inc cronos_lib dispInstall
#	@if test ! -e $(libpath)/$(libname); then \
#		echo creating symolic link to shared library $(libname); \
#		cd $(libpath) ; ln -s $(srcpath)/$(libname) $(libname); \
#	fi
	@for INC in $(incobjs); do \
		if test ! -e $(incpath)/$$INC; then \
			echo creating symbolic link to header file $$INC; \
			cd $(incpath); ln -s $(srcpath)/$$INC $$INC; \
		fi; \
	done
	@echo

cronos_lib: $(libobjs)
#	$(CPP) -shared -o $(libname) $(libobjs) -L$(thirdppath)

cronos_inc:
	@for INC in $(incobjs); do \
		if test ! -e $(incpath)/$$INC; then \
			echo creating symbolic link to header file $$INC; \
			cd $(incpath); ln -s $(srcpath)/$$INC $$INC; \
		fi; \
	done

%.o : %.cpp
	$(CPP) -c $(FLAG_CPP) $(INC_DEP) $< -o $@

%.o : %.c
	$(CPP) -c $(FLAG_CPP) $(INC_DEP) $< -o $@

%.c : $(PATH_GAMS)/apifiles/C/api/%.c
	cp $< $@

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
	-(cd $(thirdppath); make -f sobol.mk clean)
	rm -f $(libobjs) $(libname)
	-(cd $(incpath) ; rm -f $(incobjs))
	-(cd $(libpath) ; rm -f $(libname))

dispCleanInstall:
	@echo
	@(echo '***Uninstalling CRONOS library (ver.' $(version)')***')
	@echo
