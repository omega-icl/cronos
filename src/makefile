# This is the main makefile for CRONOS
# - 'make install' creates the directories 'include' and 'lib' (if not already
#   existing), compiles a shared library, and creates symbolic links to header files
#   and the compiled libraries
# - `make cleaninstall' removes the libraries, object files, and links

export version = 4.0

export rootpath       = $(shell cd .. ; pwd)
export srcpath        = $(rootpath)/src
export interfacepath  = $(rootpath)/src/interface
export incpath        = $(rootpath)/include
export binpath        = $(rootpath)/bin
export libpath        = $(rootpath)/lib

#####

all: dirs_create cronos_install interface_install

dirs_create :
	@if test ! -e $(incpath)/; then \
		echo creating directory for header files; \
		mkdir $(incpath); \
	fi
	@if test ! -e $(libpath)/; then \
		echo creating directory for compiled libraries; \
		mkdir $(libpath); \
	fi

cronos_install:
	-(cd $(srcpath); make -f cronos.mk install)

interface_install:
	-(cd $(interfacepath); make -f cronos.mk install)

#####

clean: cronos_clean interface_clean

cronos_clean:
	-(cd $(srcpath); make -f cronos.mk clean)

interface_clean:
	-(cd $(interfacepath); make -f cronos.mk clean)

#####

cleaninstall: cronos_cleandist interface_cleandist

cronos_cleandist:
	-(cd $(srcpath); make -f cronos.mk cleandist)

interface_cleandist:
	-(cd $(interfacepath); make -f cronos.mk cleandist)

