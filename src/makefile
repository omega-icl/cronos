# This is the main makefile for CRONOS
# Run `make help' for a list of available targets

export version = 4.0

export rootpath       = $(shell cd .. ; pwd)
export srcpath        = $(rootpath)/src
export interfacepath  = $(rootpath)/src/interface
export incpath        = $(rootpath)/include
export binpath        = $(rootpath)/bin
export libpath        = $(rootpath)/lib

#####

install: dirs_create cronos_install interface_install

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

uninstall: cronos_uninstall interface_uninstall

cronos_uninstall:
	-(cd $(srcpath); make -f cronos.mk uninstall)

interface_uninstall:
	-(cd $(interfacepath); make -f cronos.mk uninstall)

#####

help:
	@echo
	@echo "CRONOS Makefile (ver. $(version))"
	@echo "Available targets:"
	@echo "  install      - Create 'include' and 'lib' directories, compile a shared library,"
	@echo "                 create symbolic links to header files and compiled libraries"
	@echo "  clean        - Remove object files"
	@echo "  uninstall    - Remove libraries, object files, and symbolic links"
	@echo "  help         - Show this help message"
	@echo
