# This is the main makefile for CRONOS
# - 'make install' creates the directories 'include' and 'lib' (if not already
#   existing), compiles a shared library, and creates symbolic links to header files
#   and the compiled libraries
# - `make cleaninstall' removes the libraries, object files, and links

export version = 2.0

export rootpath = $(shell cd .. ; pwd)
export srcpath  = $(rootpath)/src
export incpath  = $(rootpath)/include
export binpath  = $(rootpath)/bin
export libpath  = $(rootpath)/lib

#####

install: dirs cronos_inst

dirs :
	@if test ! -e $(incpath)/; then \
		echo creating directory for header files; \
		mkdir $(incpath); \
	fi
	@if test ! -e $(libpath)/; then \
		echo creating directory for compiled libraries; \
		mkdir $(libpath); \
	fi
#	@if test ! -e $(binpath)/; then \
#		echo creating directory for executables; \
#		mkdir $(binpath); \
#	fi

cronos_inst:
	-(cd $(srcpath); make -f cronos.mk install)

#####

clean: cronos_clean

cronos_clean:
	-(cd $(srcpath); make -f cronos.mk clean)

#####

cleaninstall: cronos_cleandist

cronos_cleandist:
	-(cd $(srcpath); make -f cronos.mk cleandist)

