#############################################################################
#
# name: makefile
# date: 17 Sep 14
# auth: Zach Hartwig
# mail: hartwig@psfc.mit.edu
#
# desc: This file is the GNU makefile that controls the ZKBrem
#       build system. Users beware: this is no ordinary Geant4 code!
#       The build system handles a number fancy maneuvers including:
#
#       -- optional parallel build of ZKBrem with Open MPI 
#       -- generating  dictionaries for data readout into ROOT framework 
#
# dpnd: 0. The ROOT toolkit     (mandatory)
#       1. Geant4 build with Qt (optional)
#       2. Open MPI             (optional)
#
#############################################################################

# Sequential and parallel binary names
SEQ_TARGET := ZKBrem
PAR_TARGET := ZKBrem_MPI

ACROLIBS := libZKBremRoot.so

# Architecture: sequential (default), parallel (set by .PHONY 'parallel')
ARCH := sequential


##################
#  G4 Makefile  #
##################

# Set the architecture-dependent 'name' and G4TARGET variables. Note
# that the 'name' variable is required to set the name of the build
# directory in $G4INSTALL/tmp/$G4SYSTEM directory.

ifeq ($(ARCH),sequential)
	name := $(SEQ_TARGET)
else
	name := $(PAR_TARGET)
endif

G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

 .PHONY: all
all: lib/$(ACROLIBS) lib bin

# Exclude visualization for parallelized builds

ifeq ($(ARCH),parallel)
    include $(G4INSTALL)/config/binmake.gmk#.NO_VISUALIZATION
else
    include $(G4INSTALL)/config/binmake.gmk
endif	

CXX = g++

# Deal with Geant4 version parsing.  Inelegant, but should work
G4V := $(shell geant4-config --version)
G4VS := $(subst ., ,$(G4V))
G4MAJV := $(word 1,$(G4VS))
G4MINV := $(word 2,$(G4VS))
G4SUBV := $(word 3,$(G4VS))
CXXFLAGS := -DG4MAJV=$(G4MAJV) -DG4MINV=$(G4MINV) -DG4SUBV=$(G4SUBV)

# Newest version of G4 with ROOT results in massive warnings output
# since ROOT local variable 's' shadows the G4Unit 's' for
# seconds. This doesn't effect anything so suppress warning.
CXXFLAGS += $(subst -Wshadow,,$(CXXFLAGS))

CXXFLAGS += -O3 -std=c++11
CXXFLAGS += -fPIC


# Add options for debugging
#CXXFLAGS += -g
#CXXFLAGS += -fno-omit-frame-pointer
#CXXFLAGS += -O0

# turn off some clang warnings
#CXXFLAGS += -Wno-infinite-recursion -Wno-undefined-var-template

##########
#  ROOT  #
##########

# ROOT classes are presently used in ZKBrem for their immense
# utility; this requires compiling and linking against ROOT as a
# dependency. Use 'root-config' to obtain the header and library dirs
ROOTINCLUDES = -I$(shell root-config --incdir)
ROOTDISTLIBS = $(shell root-config --nonew --libs --glibs)

CPPFLAGS += $(ROOTINCLUDES)
LDLIBS +=  $(ROOTDISTLIBS) -L./lib -lZKBremRoot

ZK_ROOT_FILES = $(wildcard include/*.rhh)
ZK_ROOT_FILES += include/runMetadata.hh

lib/ZKBremDict.o : lib/ZKBremDict.cc
	@echo -e "\nBuilding $@ ...\n"
	@$(CXX) -fPIC $(CXXFLAGS) $(ROOTINCLUDES) -I. -c -o $@ $<

lib/libZKBremRoot.so : lib/ZKBremDict.o 
	@echo -e "\nBuilding the ZKBrem ROOT library ...\n"
	@$(CXX) $(ROOTDISTLIBS) -shared -o $@  $^

lib/ZKBremDict.cc : $(ZK_ROOT_FILES) include/RootLinkDef.hh
	@echo -e "Generating the ZKBrem ROOT dictionary ..."
	@rootcint -f $@ -c $^


#########################
#  Parallization / MPI  #
#########################

# Necessary flags to build ZKBrem with MPI parallelization via the
# customized acroMPImanager/messenger interface.  This conditional is
# triggered automatically via the "parallel.sh" build script.
ifeq ($(ARCH),parallel)

	# Set the Open MPI C++ compiler name
	CXX := mpic++

	# Flag for parallel build enabled by parallel.sh 
	ZK_MPI_ENABLED := 1

	# Necessary flags for parallel compilation
	CPPFLAGS += -I$(ZKBREM_MPIHOME)/include/ \
              -DZK_MPI_ENABLED
endif


.PHONY:

install:
	@echo -e "\nInstalling the ZKBrem libary in $(ZKBREM_TOPDIR)/lib/ ..."
	@cp -v $(G4WORKDIR)/tmp/$(G4SYSTEM)/ZKBrem/libZKBrem.so $(ZKBREM_TOPDIR)/lib/

libclean:
	@echo -e "\nCleaning up the ZKBrem libraries ...\n"
	@rm -f lib/*
