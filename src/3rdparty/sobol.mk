# This makefile compiles a chared library of SOBOL and it creates symbolic links
# to the library in $(libpath) and to the header files in $(incpath)

include $(srcpath)/makeoptions.mk

#####

incobjs = sobol.hpp
libobjs = sobol.o
libname = libsobol.a

#####

lib: $(libobjs)
	ar -cvq $(libname) $(libobjs)

%.o : %.cpp
	$(CPP) -c -Ofast $(FLAG_CPP) $(INC_DEP) $< -o $@

#####

clean: 
	rm -f $(libobjs) $(libname)

