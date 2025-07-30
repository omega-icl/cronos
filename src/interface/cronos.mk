# This makefile compiles the python interface and creates a symbolic link
# to the libray in $(libpath)

include $(srcpath)/makeoptions.mk

#####

libobjs = odeslvs.o ffode.o cronos.o
libname = cronos.so
libdep  = pymc.so

#####

install: dispBuild $(libname) dispInstall
	@if test ! -e $(libpath)/$(libname); then \
		echo creating symolic link to shared library $(libname); \
		cd $(libpath) ; ln -s $(interfacepath)/$(libname); \
	fi
	@if test ! -e $(libdep); then \
		ln -s $(PATH_MC)/lib/$(libdep); \
	fi
	@if test ! -e $(libpath)/$(libdep); then \
		echo creating symolic link to shared library $(libdep); \
		cd $(libpath) ; ln -s $(PATH_MC)/lib/$(libdep); \
	fi
	@echo

$(libname): $(libobjs)
	$(CPP) -shared -Wl,--export-dynamic $(libobjs) $(LIB_CRONOS) -o $(libname) 

%.o : %.cpp
	$(CPP) $(FLAG_CPP) $(FLAG_CRONOS) $(INC_CRONOS) $(INC_PYBIND11) -I$(incpath) -c $< -o $@

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
	rm -fi $(libobjs) $(libname) $(libdep)
	rm *.dot

dispClean:
	@echo
	@(echo '***Cleaning CRONOS directory (ver.' $(version)')***')
	@echo

#####

cleandist: dispCleanInstall
	rm -fi $(libobjs) $(libname) $(libdeo)
	-(cd $(libpath) ; rm -f $(libname) $(libdep))

dispCleanInstall:
	@echo
	@(echo '***Uninstalling CRONOS library (ver.' $(version)')***')
	@echo

