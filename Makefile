#SHELL = /bin/sh
CUR_DIR := $(shell pwd)

# ====================================
#        Package Configuration
# ====================================
# Enable packages by setting to 1 on command line:
#   make PYTHON=1 AENET=1 LAMMPS=1
#
# Or set compiler:
#   make FC=gfortran PYTHON=1
#   make COMPILER=gfortran PYTHON=1
#
# Debug mode:
#   make DEBUG=1 PYTHON=1
#
# Build shared library:
#   make lib PYTHON=1

PYTHON ?= 0
AENET ?= 0
LAMMPS ?= 0
DEBUG ?= 0
COMPILER ?= mpif90

# External library paths (override on command line if needed)
LAMMPS_LIB_PATH ?= /usr/local/lib
AENET_LIB_PATH ?= $(LIB)

# ====================================
#        Compiler Options
# ====================================
FC := mpif90
AR := ar
CC := mpicc

# Intel Fortran flags
OPTIMIZE_FLAGS_IFORT := -O3
OPTIMIZE_FLAGS_IFORT += -xHost
OPTIMIZE_FLAGS_IFORT += -no-prec-div
OPTIMIZE_FLAGS_IFORT += -no-wrap-margin
OPTIMIZE_FLAGS_IFORT += -fpp
OPTIMIZE_FLAGS_IFORT += -traceback 

DEBUG_FLAGS_IFORT := -check all -traceback -g -fpe0 -O0 -fp-stack-check -debug all -ftrapuv -fpp -no-wrap-margin

# GFortran flags
OPTIMIZE_FLAGS_GFORT := -O3 -cpp -g 
OPTIMIZE_FLAGS_GFORT += -fbacktrace -fcheck=bounds -ffree-line-length-512
OPTIMIZE_FLAGS_GFORT += -finit-real=zero -finit-integer=0

DEBUG_FLAGS_GFORT := -cpp -fbacktrace -fcheck=all -g -ffree-line-length-512 -ffpe-trap=overflow,invalid,zero -Waliasing -Wsurprising

# Library build flags
LIBRARY_FLAGS := -shared -fpic

# Base linker flags
LDFLAGS := -lfftw3

# Base package flags (always enabled)
PACKAGE_FLAGS := -DMPIPARALLEL

# ====================================
#        Directory List
# ====================================
SRC := $(CUR_DIR)/src
LIB := $(CUR_DIR)/lib
MOD := $(CUR_DIR)/mods
OBJ := $(CUR_DIR)/objects
PYTHON_DIR := $(CUR_DIR)/python

# Define file extensions
.SUFFIXES:
.SUFFIXES: .f .f90 .o .mod 

# ====================================
#        Source Files
# ====================================
SRC_MAIN := $(SRC)/Common.f90\
        		$(SRC)/ErrorChecking.f90\
        		$(SRC)/C_To_Fotran.f90\
        		$(SRC)/fftw3_interface.f90\
        		$(SRC)/Common_BoxData.f90\
        		$(SRC)/Common_TrajData.f90\
        		$(SRC)/Common_Analysis.f90\
        		$(SRC)/Common_ECalc.f90\
        		$(SRC)/Common_MolDef.f90\
        		$(SRC)/Common_Sampling.f90\
        		$(SRC)/Common_MCMoves.f90\
        		$(SRC)/Common_NeighList.f90\
        		$(SRC)/Debug.f90\
        		$(SRC)/Data_Graph.f90\
        		$(SRC)/Data_Queue.f90\
         		$(SRC)/Constrain_MultiAtomDistCrit.f90\
         		$(SRC)/Constrain_MolTotal.f90\
         		$(SRC)/Constrain_DistCriteria.f90\
         		$(SRC)/Constrain_EnergyCeiling.f90\
         		$(SRC)/Constrain_EnergyFloor.f90\
         		$(SRC)/Constrain_FreezeType.f90\
	        	$(SRC)/SearchSort.f90\
        		$(SRC)/Sampling_AcceptAll.f90\
        		$(SRC)/Sampling_AcceptNone.f90\
        		$(SRC)/Sampling_Metropolis.f90\
        		$(SRC)/Sampling_MinMetrop.f90\
        		$(SRC)/Sampling_Nested.f90\
        		$(SRC)/Sampling_Umbrella.f90\
        		$(SRC)/Sampling_UmbrellaWHAM.f90\
        		$(SRC)/Math_CoordinateFunctions.f90\
        		$(SRC)/Move_MC_AVBMC.f90\
        		$(SRC)/Move_MC_AnisoVol.f90\
        		$(SRC)/Move_MC_AtomExchange.f90\
        		$(SRC)/Move_MC_AtomTranslation.f90\
        		$(SRC)/Move_MC_ParticleExchange.f90\
        		$(SRC)/Move_MC_BasicSwap.f90\
        		$(SRC)/Move_MC_CBMC.f90\
        		$(SRC)/Move_MC_Delete.f90\
        		$(SRC)/Move_MC_MolTranslation.f90\
        		$(SRC)/Move_MC_IsoVol.f90\
        		$(SRC)/Move_MC_PlaneRotate.f90\
        		$(SRC)/Move_MC_PlaneTranslate.f90\
        		$(SRC)/Move_MC_PlaneAtomTranslate.f90\
        		$(SRC)/Move_MC_UBSwap.f90\
        		$(SRC)/Move_MC_VolExchange.f90\
        		$(SRC)/ExeptionHandling.f90\
        		$(SRC)/MolSearch.f90\
        		$(SRC)/Analysis_AngleDistribution.f90\
        		$(SRC)/Analysis_BondDistribution.f90\
        		$(SRC)/Analysis_TorsionDistribution.f90\
        		$(SRC)/Analysis_TotalSize.f90\
        		$(SRC)/Analysis_BlockAverage.f90\
        		$(SRC)/Analysis_DensityOfStates.f90\
        		$(SRC)/Analysis_ClusterSize.f90\
        		$(SRC)/Analysis_DistPair.f90\
        		$(SRC)/Analysis_MolFraction.f90\
        		$(SRC)/Analysis_RDF.f90\
        		$(SRC)/Analysis_ThermoAverage.f90\
        		$(SRC)/Analysis_ThermoIntegration.f90\
        		$(SRC)/Box_Presets.f90\
        		$(SRC)/Box_CubicBox.f90\
        		$(SRC)/Box_OrthoBox.f90\
        		$(SRC)/FF_AENet.f90\
        		$(SRC)/FF_EasyPair_Cut.f90\
        		$(SRC)/FF_Lammps.f90\
        		$(SRC)/FF_Einstein.f90\
        		$(SRC)/FF_HardSphere.f90\
        		$(SRC)/FF_Hybrid.f90\
        		$(SRC)/FF_EP_LJ_Ele_Cut.f90\
        		$(SRC)/FF_EP_LJ_Cut.f90\
        		$(SRC)/FF_EP_LJ_CutShift.f90\
        		$(SRC)/FF_EP_Pedone_Cut.f90\
        		$(SRC)/FF_EP_TosiFumi_Cut.f90\
        		$(SRC)/FF_Pedone.f90\
        		$(SRC)/FF_LJ_Cut.f90\
        		$(SRC)/FF_Tersoff.f90\
        		$(SRC)/FF_StilWeb.f90\
        		$(SRC)/FF_ThermoInt.f90\
        		$(SRC)/FF_Ewald.f90\
        		$(SRC)/FF_LJ_Ewald.f90\
        		$(SRC)/FF_P3M.f90\
        		$(SRC)/FF_FMM.f90\
        		$(SRC)/SphericalHarmonics.f90\
        		$(SRC)/FF_EAM.f90\
        		$(SRC)/Intra_AngleHarmonic.f90\
        		$(SRC)/Intra_AngleRidgid.f90\
        		$(SRC)/Intra_BondRidgid.f90\
        		$(SRC)/Intra_BondHarmonic.f90\
        		$(SRC)/Intra_TorsionTrappe.f90\
        		$(SRC)/Intra_TorsionHarmonic.f90\
        		$(SRC)/Intra_TorsionRidgid.f90\
        		$(SRC)/Intra_Misc1_5_Pair.f90\
        		$(SRC)/MolCon_LinearCBMC.f90\
        		$(SRC)/MolCon_BranchCBMC.f90\
        		$(SRC)/MolCon_RidgidRegrowth.f90\
        		$(SRC)/MolCon_SimpleRegrowth.f90\
 	        	$(SRC)/Script_AnalysisType.f90\
 	        	$(SRC)/Script_AngleType.f90\
 	        	$(SRC)/Script_BondType.f90\
 	        	$(SRC)/Script_MiscType.f90\
 	        	$(SRC)/Script_TorsionType.f90\
 	        	$(SRC)/Script_Constraint.f90\
 	        	$(SRC)/Script_Forcefield.f90\
 	        	$(SRC)/Script_FieldType.f90\
 	        	$(SRC)/Script_LoadCoords.f90\
	        	$(SRC)/Script_MCMoves.f90\
	        	$(SRC)/Script_Main.f90\
 	        	$(SRC)/Script_NeighType.f90\
 	        	$(SRC)/Script_RegrowType.f90\
	        	$(SRC)/Script_Sampling.f90\
	        	$(SRC)/Script_SimBoxes.f90\
	        	$(SRC)/Script_Initialize.f90\
 	        	$(SRC)/Script_TrajType.f90\
	        	$(SRC)/Output_DumpCoords.f90\
	        	$(SRC)/Traj_POSCAR.f90\
	        	$(SRC)/Traj_LAMMPSDump.f90\
	        	$(SRC)/Traj_XYZFormat.f90\
	        	$(SRC)/Traj_XSF.f90\
	        	$(SRC)/Input_Format.f90\
 	        	$(SRC)/Neigh_CellRSqList.f90\
 	        	$(SRC)/Neigh_CellList.f90\
 	        	$(SRC)/Neigh_RSqList.f90\
        		$(SRC)/VariablePrecision.f90\
        		$(SRC)/Sim_Minimize.f90\
        		$(SRC)/Sim_MonteCarlo.f90\
        		$(SRC)/Sim_GeneticAlgor.f90\
        		$(SRC)/Sim_LibControl.f90\
        		$(SRC)/Units.f90

SRC_TEMPLATE := $(SRC)/Template_Master.f90\
        		$(SRC)/RandomNew.f90\
        		$(SRC)/Template_AcceptRule.f90\
				$(SRC)/Template_Analysis.f90\
				$(SRC)/Template_BondFF.f90\
				$(SRC)/Template_AngleFF.f90\
				$(SRC)/Template_TorsionFF.f90\
				$(SRC)/Template_MiscIntra.f90\
				$(SRC)/Template_Constraint.f90\
				$(SRC)/Template_Forcefield.f90\
				$(SRC)/Template_SimBox.f90\
        $(SRC)/Box_SimpleBox.f90\
				$(SRC)/Box_Ultility.f90\
				$(SRC)/Template_IntraFF.f90\
				$(SRC)/Template_NeighList.f90\
				$(SRC)/Template_Trajectory.f90\
				$(SRC)/Template_MultiBoxMove.f90\
				$(SRC)/Template_MoveClass.f90\
				$(SRC)/Template_MolConstructor.f90

# ====================================
#        Python Package Sources
# ====================================
# These modules are always compiled (have stubs when EMBPYTHON not defined)
# because Script files reference them unconditionally
SRC_PYTHON_ALWAYS := $(PYTHON_DIR)/Traj_Python.f90\
                     $(PYTHON_DIR)/Constrain_Python.f90

# These require full Python embedding support
SRC_PYTHON := $(PYTHON_DIR)/forpy_mod.f90\
              $(PYTHON_DIR)/Python_CommonTypes.f90\
              $(PYTHON_DIR)/Analysis_Python.f90\
              $(PYTHON_DIR)/Move_Python.f90\
              $(PYTHON_DIR)/Energy_Python.f90\
              $(PYTHON_DIR)/Sim_Python.f90

# ====================================
#        Object Files
# ====================================
OBJ_LIBRARY = $(shell find $(LIB) -name '*.a' 2>/dev/null)

# Include any additional makefiles (for backward compatibility)
-include *.Makefile 

# Always include Traj_Python (has stubs for non-Python builds)
SRC_COMPLETE := $(SRC_TEMPLATE) $(SRC_MAIN) $(SRC_PYTHON_ALWAYS)

OBJ_MAIN := $(patsubst $(SRC)/%.f90, $(OBJ)/%.o, $(SRC_MAIN))
OBJ_MAIN := $(patsubst $(PYTHON_DIR)/%.f90, $(OBJ)/%.o, $(OBJ_MAIN))
OBJ_TEMPLATE := $(patsubst $(SRC)/%.f90, $(OBJ)/%.o, $(SRC_TEMPLATE))
OBJ_PYTHON := $(patsubst $(PYTHON_DIR)/%.f90, $(OBJ)/%.o, $(SRC_PYTHON))
OBJ_PYTHON_ALWAYS := $(patsubst $(PYTHON_DIR)/%.f90, $(OBJ)/%.o, $(SRC_PYTHON_ALWAYS))
OBJ_COMPLETE := $(OBJ_TEMPLATE) $(OBJ_MAIN) $(OBJ_PYTHON_ALWAYS) 

# ====================================
#        Package Configuration Logic
# ====================================

# Select compiler flags based on COMPILER variable
ifeq ($(COMPILER),gfortran)
  ifeq ($(DEBUG),1)
    BASE_COMPFLAGS := $(DEBUG_FLAGS_GFORT)
  else
    BASE_COMPFLAGS := $(OPTIMIZE_FLAGS_GFORT)
  endif
else
  # Default to Intel
  ifeq ($(DEBUG),1)
    BASE_COMPFLAGS := $(DEBUG_FLAGS_IFORT)
  else
    BASE_COMPFLAGS := $(OPTIMIZE_FLAGS_IFORT)
  endif
endif

# Python package configuration
PYTHON_INCLUDE :=
PYTHON_LDFLAGS :=
PYTHON_LIBS :=

ifeq ($(PYTHON),1)
  PACKAGE_FLAGS += -DEMBPYTHON
  
  # Python configuration (auto-detected)
  # Python 3.8+ requires --embed flag for embedding Python
  PYTHON_CONFIG := $(shell which python3-config 2>/dev/null || which python-config 2>/dev/null)
  ifneq ($(PYTHON_CONFIG),)
    PYTHON_INCLUDE := $(shell $(PYTHON_CONFIG) --includes 2>/dev/null)
    # Try --embed first (Python 3.8+), fallback to without
    PYTHON_LDFLAGS := $(shell $(PYTHON_CONFIG) --ldflags --embed 2>/dev/null || $(PYTHON_CONFIG) --ldflags 2>/dev/null)
    PYTHON_LIBS := $(shell $(PYTHON_CONFIG) --libs --embed 2>/dev/null || $(PYTHON_CONFIG) --libs 2>/dev/null)
  else
    # Fallback to sysconfig-based detection
    PYTHON_INCLUDE := -I$(shell python3 -c "import sysconfig; print(sysconfig.get_path('include'))" 2>/dev/null)
    PYTHON_LIBDIR := $(shell python3 -c "import sysconfig; print(sysconfig.get_config_var('LIBDIR'))" 2>/dev/null)
    PYTHON_VERSION := $(shell python3 -c "import sysconfig; print(sysconfig.get_config_var('LDVERSION') or sysconfig.get_config_var('VERSION'))" 2>/dev/null)
    PYTHON_LDFLAGS := -L$(PYTHON_LIBDIR)
    PYTHON_LIBS := -lpython$(PYTHON_VERSION)
  endif
  
  LDFLAGS += $(PYTHON_LDFLAGS) $(PYTHON_LIBS)
  
  # Add Python sources
  SRC_COMPLETE := $(SRC_PYTHON) $(SRC_COMPLETE)
  OBJ_COMPLETE := $(OBJ_PYTHON) $(OBJ_COMPLETE)
endif

# AENet package configuration
ifeq ($(AENET),1)
  PACKAGE_FLAGS += -DAENET
  LDFLAGS += -llapack -lblas
  ifneq ($(wildcard $(AENET_LIB_PATH)/libaenet.a),)
    OBJ_LIBRARY += $(AENET_LIB_PATH)/libaenet.a
  endif
endif

# LAMMPS package configuration
ifeq ($(LAMMPS),1)
  PACKAGE_FLAGS += -DLAMMPS
  LDFLAGS += -L$(LAMMPS_LIB_PATH) -llammps -lstdc++
endif

# Build final COMPFLAGS with all package flags and Python includes
COMPFLAGS := $(BASE_COMPFLAGS) $(PACKAGE_FLAGS) $(PYTHON_INCLUDE)

# ====================================
#        Build Targets
# ====================================

# Default target
default: startUP classyMC modout finale

# Shared library target
lib: COMPFLAGS += $(LIBRARY_FLAGS)
lib: startUP libclassymc.so modout finale

# Debug target (can also use DEBUG=1)
debug: DEBUG=1
debug: COMPFLAGS := $(DEBUG_FLAGS_IFORT) $(PACKAGE_FLAGS)
debug: startUP classyMC modout finale

# GFortran target (can also use COMPILER=gfortran)
gfortran: COMPILER=gfortran
gfortran: COMPFLAGS := $(OPTIMIZE_FLAGS_GFORT) $(PACKAGE_FLAGS)
gfortran: startUP classyMC modout finale

# Legacy targets for backward compatibility
aenet: AENET=1
aenet: COMPFLAGS := $(BASE_COMPFLAGS) $(PACKAGE_FLAGS) -DAENET
aenet: LDFLAGS += -llapack -lblas
aenet: startUP classyMC modout finale

lammps: LAMMPS=1
lammps: COMPFLAGS := $(BASE_COMPFLAGS) $(PACKAGE_FLAGS) -DLAMMPS
lammps: LDFLAGS += -L$(LAMMPS_LIB_PATH) -llammps -lstdc++
lammps: startUP classyMC modout finale

# Clean targets
clean: removeObjects removeExec finale

# Help target
help:
	@echo "=============================================================="
	@echo "                    ClassyMC Build System"
	@echo "=============================================================="
	@echo ""
	@echo "Usage: make [target] [OPTIONS]"
	@echo ""
	@echo "Targets:"
	@echo "  default     Build ClassyMC executable (default)"
	@echo "  lib         Build shared library (libclassymc.so)"
	@echo "  debug       Build with debug flags"
	@echo "  gfortran    Build with GFortran compiler"
	@echo "  clean       Remove all build artifacts"
	@echo "  help        Show this help message"
	@echo ""
	@echo "Package Options (set to 1 to enable):"
	@echo "  PYTHON=1    Enable embedded Python support"
	@echo "              Allows Python scripts for moves, analysis, trajectories"
	@echo ""
	@echo "  AENET=1     Enable AENet neural network potentials"
	@echo "              Requires LAPACK/BLAS libraries"
	@echo ""
	@echo "  LAMMPS=1    Enable LAMMPS as a force field backend"
	@echo "              Set LAMMPS_LIB_PATH if not in /usr/local/lib"
	@echo ""
	@echo "Build Options:"
	@echo "  DEBUG=1          Enable debug mode with detailed checks"
	@echo "  COMPILER=ifort   Use Intel Fortran (default)"
	@echo "  COMPILER=gfortran Use GNU Fortran"
	@echo ""
	@echo "Library Paths (override defaults):"
	@echo "  LAMMPS_LIB_PATH=/path/to/lammps/lib"
	@echo "  AENET_LIB_PATH=/path/to/aenet/lib"
	@echo ""
	@echo "Examples:"
	@echo "  make                           # Basic build"
	@echo "  make PYTHON=1                  # Build with Python support"
	@echo "  make PYTHON=1 AENET=1          # Python + AENet"
	@echo "  make DEBUG=1 PYTHON=1          # Debug build with Python"
	@echo "  make COMPILER=gfortran PYTHON=1  # GFortran with Python"
	@echo "  make lib PYTHON=1              # Shared library with Python"
	@echo "  make LAMMPS=1 LAMMPS_LIB_PATH=/opt/lammps/lib"
	@echo ""
	@echo "=============================================================="

# ====================================
#        Compile Rules
# ====================================

%.o: %.mod

.f90.o:
	@echo "Compiling $<"
	@$(FC) $(COMPFLAGS) $(MODFLAGS) -c -o $@ $<

.f90.mod:
	@echo "Compiling $<"
	@$(FC) $(COMPFLAGS) $(MODFLAGS) -c -o $@ $<

$(OBJ)/%.o: $(SRC)/%.f90
	@echo "Compiling $<"
	@$(FC) $(COMPFLAGS) $(MODFLAGS) -c $< -o $@ 

$(OBJ)/%.o: $(PYTHON_DIR)/%.f90
	@echo "Compiling $<"
	@$(FC) $(COMPFLAGS) $(MODFLAGS) -c $< -o $@ 

$(OBJ)/%.o: $(TEMPLATE)/%.f90
	@echo "Compiling $<"
	@$(FC) $(COMPFLAGS) $(MODFLAGS) -c $< -o $@ 

# ====================================
#        Link Targets
# ====================================

libclassymc.so: $(OBJ_COMPLETE) $(OBJ_LIBRARY)
	@echo "============================================="
	@echo "    Linking Shared Library"
	@echo "============================================="
	@$(FC) $(COMPFLAGS) $(LIBRARY_FLAGS) $(MODFLAGS) $^ -o $@ $(LDFLAGS)

classyMC: $(OBJ_COMPLETE) $(SRC)/Main.f90 $(OBJ_LIBRARY)
	@echo "============================================="
	@echo "    Linking ClassyMC Executable"
	@echo "============================================="
	@echo "    Packages enabled:"
ifeq ($(PYTHON),1)
	@echo "      - Python (embedded interpreter)"
endif
ifeq ($(AENET),1)
	@echo "      - AENet (neural network potentials)"
endif
ifeq ($(LAMMPS),1)
	@echo "      - LAMMPS (force field backend)"
endif
	@echo "============================================="
	@$(FC) $(COMPFLAGS) $(MODFLAGS) $^ -o $@ $(LDFLAGS)

# ====================================
#        Utility Targets
# ====================================

startUP:
	@echo "=================================================================="
	@echo "                     ClassyMC Build Starting"
	@echo "=================================================================="
	@echo "  Directory: $(CUR_DIR)"
	@echo "  Compiler:  $(FC)"
	@echo "  Flags:     $(COMPFLAGS)"
ifeq ($(PYTHON),1)
	@echo "  Python:    ENABLED"
	@echo "    Include: $(PYTHON_INCLUDE)"
	@echo "    Libs:    $(PYTHON_LIBS)"
endif
ifeq ($(AENET),1)
	@echo "  AENet:     ENABLED"
endif
ifeq ($(LAMMPS),1)
	@echo "  LAMMPS:    ENABLED ($(LAMMPS_LIB_PATH))"
endif
	@echo "=================================================================="
	@mv $(MOD)/*.mod $(CUR_DIR)/ 2>/dev/null || true

modout:
	@mv $(CUR_DIR)/*.mod $(MOD)/ 2>/dev/null || true

finale:
	@echo ""
	@echo "======================== Build Complete ========================="
	@echo "=================================================================="

removeObjects:
	@echo "============================================="
	@echo "         Cleaning Build Artifacts"
	@echo "============================================="
	@rm -f ./*.o ./*.mod
	@rm -f $(SRC)/*.o $(SRC)/*.mod
	@rm -f $(MOD)/*.o $(MOD)/*.mod
	@rm -f $(SRC)/*/*.o $(SRC)/*/*.mod
	@rm -f $(OBJ)/*.o

removeExec:
	@rm -f $(CUR_DIR)/Python.Makefile
	@rm -f $(CUR_DIR)/libclassymc.so
	@rm -f $(CUR_DIR)/classyMC
	@rm -f $(CUR_DIR)/classyMCAENet
	@rm -f $(CUR_DIR)/classyMCLAMMPS
	@rm -f $(CUR_DIR)/classyMC_debug
	@rm -f $(CUR_DIR)/classyMC.exe

# ====================================
#        Dependencies
# ====================================
$(OBJ)/Common.o: $(OBJ)/VariablePrecision.o
$(OBJ)/Common_BoxData.o: $(OBJ)/Box_SimpleBox.o 
$(OBJ)/Common_Analysis.o: $(OBJ)/Template_Analysis.o
$(OBJ)/Common_ECalc.o: $(OBJ)/Template_Forcefield.o $(OBJ)/Common.o
$(OBJ)/Common_Sampling.o: $(OBJ)/Template_AcceptRule.o $(OBJ)/Sampling_Metropolis.o
$(OBJ)/Common_MolDef.o: $(OBJ)/Template_MolConstructor.o $(OBJ)/Template_BondFF.o $(OBJ)/Template_TorsionFF.o $(OBJ)/Template_AngleFF.o  $(OBJ)/Template_MiscIntra.o $(OBJ)/Data_Graph.o
${OBJ}/Common_TrajData.o: $(OBJ)/Template_Trajectory.o

$(OBJ)/Template_Constraint.o: $(OBJ)/Template_SimBox.o $(OBJ)/Template_Master.o
$(OBJ)/Template_AcceptRule.o: $(OBJ)/Common.o $(OBJ)/Input_Format.o $(OBJ)/Template_SimBox.o $(OBJ)/Template_Master.o
$(OBJ)/Template_Anaylsis.o: $(OBJ)/Box_SimpleBox.o $(OBJ)/Template_Master.o
$(OBJ)/Template_SimBox.o: $(OBJ)/Common.o ${OBJ}/Input_Format.o $(OBJ)/Template_NeighList.o $(OBJ)/Template_Master.o
$(OBJ)/Template_Master.o: $(OBJ)/VariablePrecision.o
$(OBJ)/Template_MultiBoxMove.o: $(OBJ)/Template_MoveClass.o $(OBJ)/Template_Master.o $(OBJ)/Box_Ultility.o
$(OBJ)/Template_MoveClass.o: $(OBJ)/Common.o $(OBJ)/Common_BoxData.o ${OBJ}/Box_SimpleBox.o $(OBJ)/Template_Master.o $(OBJ)/Box_Ultility.o $(OBJ)/RandomNew.o
$(OBJ)/Template_Forcefield.o: $(OBJ)/Common.o $(OBJ)/Common_MolDef.o $(OBJ)/Template_SimBox.o $(OBJ)/Template_Master.o
$(OBJ)/Template_NeighList.o: $(OBJ)/SearchSort.o $(OBJ)/Template_Master.o $(OBJ)/Input_Format.o
$(OBJ)/Template_MolConstructor.o: $(OBJ)/Template_SimBox.o $(OBJ)/Template_Master.o
$(OBJ)/Template_AngleFF.o: $(OBJ)/Template_IntraFF.o $(OBJ)/Template_Master.o
$(OBJ)/Template_BondFF.o: $(OBJ)/Template_IntraFF.o $(OBJ)/Template_Master.o $(OBJ)/RandomNew.o
$(OBJ)/Template_TorsionFF.o: $(OBJ)/Template_IntraFF.o $(OBJ)/Template_Master.o

$(OBJ)/Analysis_ThermoIntegration.o: $(OBJ)/FF_ThermoInt.o
$(OBJ)/Box_SimpleBox.o: $(OBJ)/Common.o $(OBJ)/Template_NeighList.o $(OBJ)/Input_Format.o $(OBJ)/Common_ECalc.o $(OBJ)/Template_SimBox.o $(OBJ)/Template_Constraint.o $(OBJ)/Units.o $(OBJ)/Common_NeighList.o $(OBJ)/ErrorChecking.o
$(OBJ)/Box_CubicBox.o: $(OBJ)/Box_SimpleBox.o
$(OBJ)/Box_OrthoBox.o: $(OBJ)/Box_SimpleBox.o
$(OBJ)/Box_Ultility.o: $(OBJ)/Box_SimpleBox.o
$(OBJ)/Box_Presets.o: $(OBJ)/Box_OrthoBox.o $(OBJ)/Box_CubicBox.o

$(OBJ)/Move_MC_AVBMC.o: $(OBJ)/Common.o $(OBJ)/Box_Ultility.o
$(OBJ)/Move_MC_CBMC.o: $(OBJ)/Common.o $(OBJ)/Box_Ultility.o $(OBJ)/MolCon_LinearCBMC.o $(OBJ)/MolCon_BranchCBMC.o  
$(OBJ)/Move_MC_AtomTranslation.o: $(OBJ)/Common.o $(OBJ)/Common_BoxData.o $(OBJ)/Box_SimpleBox.o $(OBJ)/RandomNew.o $(OBJ)/Template_MoveClass.o $(OBJ)/Template_Constraint.o $(OBJ)/Box_Ultility.o $(OBJ)/Common_Sampling.o $(OBJ)/Move_MC_MolTranslation.o
$(OBJ)/Move_MC_IsoVol.o: $(OBJ)/Common.o $(OBJ)/Common_BoxData.o $(OBJ)/Box_CubicBox.o $(OBJ)/Box_OrthoBox.o $(OBJ)/RandomNew.o $(OBJ)/Template_MoveClass.o $(OBJ)/Template_Constraint.o $(OBJ)/Box_Ultility.o
$(OBJ)/Move_MC_AnisoVol.o: $(OBJ)/Common.o $(OBJ)/Common_BoxData.o $(OBJ)/Box_CubicBox.o $(OBJ)/Box_OrthoBox.o $(OBJ)/RandomNew.o $(OBJ)/Template_MoveClass.o $(OBJ)/Template_Constraint.o $(OBJ)/Box_Ultility.o
$(OBJ)/Move_MC_AtomExchange.o: $(OBJ)/Common.o $(OBJ)/Common_BoxData.o $(OBJ)/Box_SimpleBox.o $(OBJ)/RandomNew.o $(OBJ)/Template_MoveClass.o $(OBJ)/Box_Ultility.o
$(OBJ)/Move_MC_ThermoLambda.o: $(OBJ)/FF_ThermoInt.o $(OBJ)/Analysis_ThermoIntegration.o 
$(OBJ)/Move_GA_AtomExchange.o: $(OBJ)/Common.o $(OBJ)/Common_BoxData.o $(OBJ)/Box_Ultility.o
$(OBJ)/Move_MC_MolTranslation.o: $(OBJ)/Common_Sampling.o
$(OBJ)/MolCon_SimpleRegrowth.o: $(OBJ)/Template_MolConstructor.o
$(OBJ)/MolCon_LinearCBMC.o: $(OBJ)/Template_MolConstructor.o $(OBJ)/MolSearch.o
$(OBJ)/MolCon_BranchCBMC.o: $(OBJ)/Template_MolConstructor.o $(OBJ)/MolSearch.o $(OBJ)/Data_Queue.o $(OBJ)/Math_CoordinateFunctions.o
$(OBJ)/Data_Queue.o: $(OBJ)/VariablePrecision.o

$(OBJ)/Script_Main.o: $(OBJ)/Units.o $(OBJ)/Common_BoxData.o $(OBJ)/Script_Forcefield.o $(OBJ)/Box_CubicBox.o $(OBJ)/Script_SimBoxes.o $(OBJ)/Script_Sampling.o $(OBJ)/Script_MCMoves.o $(OBJ)/Script_Initialize.o $(OBJ)/Script_NeighType.o $(OBJ)/Script_TrajType.o $(OBJ)/Sim_MonteCarlo.o $(OBJ)/Sim_Minimize.o

$(OBJ)/Script_Forcefield.o: ${OBJ}/Input_Format.o ${OBJ}/Template_Forcefield.o ${OBJ}/Move_MC_AtomTranslation.o ${OBJ}/Units.o $(OBJ)/Script_FieldType.o $(OBJ)/Script_BondType.o $(OBJ)/Script_AngleType.o $(OBJ)/Script_RegrowType.o 
$(OBJ)/Script_LoadCoords.o: ${OBJ}/Script_SimBoxes.o
$(OBJ)/Script_FieldType.o: ${OBJ}/Input_Format.o ${OBJ}/Template_Forcefield.o ${OBJ}/FF_LJ_Cut.o ${OBJ}/Move_MC_AtomTranslation.o $(OBJ)/Common_ECalc.o ${OBJ}/FF_Lammps.o
$(OBJ)/Script_TrajType.o: ${OBJ}/Common_TrajData.o ${OBJ}/Template_Trajectory.o ${OBJ}/Traj_XSF.o ${OBJ}/Traj_XYZFormat.o $(OBJ)/Traj_LAMMPSDump.o $(OBJ)/Traj_POSCAR.o $(OBJ)/Traj_Python.o
$(OBJ)/Script_NeighType.o: ${OBJ}/Neigh_RSqList.o $(OBJ)/Neigh_CellRSqList.o $(OBJ)/Neigh_CellList.o $(OBJ)/Common_BoxData.o

$(OBJ)/RandomNew.o: $(OBJ)/Common.o $(OBJ)/Units.o

$(OBJ)/Neigh_RSqList.o: $(OBJ)/Common_BoxData.o $(OBJ)/Template_NeighList.o $(OBJ)/Common_NeighList.o
$(OBJ)/Neigh_CellRSqList.o: $(OBJ)/Neigh_RSqList.o $(OBJ)/Common_BoxData.o
$(OBJ)/Neigh_CellList.o: $(OBJ)/Common_BoxData.o $(OBJ)/Template_NeighList.o
$(OBJ)/Sampling_Umbrella.o: $(OBJ)/Sampling_UmbrellaWHAM.o
$(OBJ)/Sampling_Metropolis.o: $(OBJ)/RandomNew.o

$(OBJ)/Main.o: $(OBJ)/Sim_MonteCarlo.o $(OBJ)/Sim_Minimize.o
$(OBJ)/Sim_Library.o: $(OBJ)/Script_Main.o

$(OBJ)/FF_EP_LJ_Cut.o: $(OBJ)/FF_EasyPair_Cut.o
$(OBJ)/FF_EP_Pedone_Cut.o: $(OBJ)/FF_EasyPair_Cut.o
$(OBJ)/FF_EP_TosiFumi_Cut.o: $(OBJ)/FF_EasyPair_Cut.o
$(OBJ)/FF_EP_LJ_CutShift.o: $(OBJ)/FF_EasyPair_Cut.o
$(OBJ)/FF_Lammps.o: $(OBJ)/FF_EasyPair_Cut.o $(OBJ)/Box_CubicBox.o $(OBJ)/Box_OrthoBox.o
$(OBJ)/FF_Ewald.o: $(OBJ)/Template_Forcefield.o $(OBJ)/Units.o $(OBJ)/Box_OrthoBox.o $(OBJ)/Box_CubicBox.o
$(OBJ)/FF_LJ_Ewald.o: $(OBJ)/Template_Forcefield.o $(OBJ)/Units.o $(OBJ)/Box_OrthoBox.o $(OBJ)/Box_CubicBox.o
$(OBJ)/FF_P3M.o: $(OBJ)/FF_EasyPair_Cut.o $(OBJ)/Template_Forcefield.o $(OBJ)/Units.o $(OBJ)/Box_OrthoBox.o $(OBJ)/Box_CubicBox.o $(OBJ)/fftw3_interface.o
$(OBJ)/FF_EAM.o: $(OBJ)/Template_Forcefield.o $(OBJ)/Units.o
$(OBJ)/SphericalHarmonics.o: $(OBJ)/VariablePrecision.o
$(OBJ)/FF_FMM.o: $(OBJ)/Template_Forcefield.o $(OBJ)/Units.o $(OBJ)/Box_OrthoBox.o $(OBJ)/Box_CubicBox.o $(OBJ)/SphericalHarmonics.o
$(OBJ)/Script_FieldType.o: $(OBJ)/FF_Ewald.o $(OBJ)/FF_LJ_Ewald.o $(OBJ)/FF_P3M.o $(OBJ)/FF_FMM.o $(OBJ)/FF_EAM.o
$(OBJ)/Move_MC_PlaneAtomTranslate.o: $(OBJ)/Move_MC_PlaneTranslate.o

$(OBJ)/Sim_MonteCarlo.o: $(OBJ)/Common.o $(OBJ)/Units.o $(OBJ)/Move_MC_AtomTranslation.o $(OBJ)/RandomNew.o $(OBJ)/Common_TrajData.o $(OBJ)/Output_DumpCoords.o $(OBJ)/Common_Analysis.o $(OBJ)/Common_MCMoves.o $(OBJ)/Template_MultiBoxMove.o

# ====================================
#        Python Package Dependencies
# ====================================
# These modules are always compiled (have stubs), base dependencies only
$(OBJ)/Traj_Python.o: $(OBJ)/Template_Trajectory.o $(OBJ)/Input_Format.o
$(OBJ)/Constrain_Python.o: $(OBJ)/Template_Constraint.o $(OBJ)/Input_Format.o $(OBJ)/Box_SimpleBox.o
$(OBJ)/Script_Constraint.o: $(OBJ)/Constrain_Python.o

# Full Python embedding dependencies (only when PYTHON=1)
ifeq ($(PYTHON),1)
$(OBJ)/Python_CommonTypes.o: $(OBJ)/forpy_mod.o $(OBJ)/Common_BoxData.o $(OBJ)/Box_CubicBox.o $(OBJ)/Box_OrthoBox.o
$(OBJ)/Sim_Python.o: $(OBJ)/forpy_mod.o $(OBJ)/Sim_MonteCarlo.o
$(OBJ)/Analysis_Python.o: $(OBJ)/forpy_mod.o $(OBJ)/Template_Analysis.o $(OBJ)/Input_Format.o $(OBJ)/Common_Analysis.o $(OBJ)/Common_BoxData.o $(OBJ)/Python_CommonTypes.o
$(OBJ)/Move_Python.o: $(OBJ)/forpy_mod.o $(OBJ)/Template_MoveClass.o $(OBJ)/Input_Format.o $(OBJ)/Common_MCMoves.o $(OBJ)/Common_BoxData.o $(OBJ)/Python_CommonTypes.o $(OBJ)/Common_Sampling.o
# Additional dependencies when Python is enabled
$(OBJ)/Traj_Python.o: $(OBJ)/forpy_mod.o $(OBJ)/Common_BoxData.o $(OBJ)/Box_CubicBox.o $(OBJ)/Box_OrthoBox.o
$(OBJ)/Constrain_Python.o: $(OBJ)/forpy_mod.o $(OBJ)/Common_BoxData.o $(OBJ)/Box_CubicBox.o $(OBJ)/Box_OrthoBox.o $(OBJ)/Python_CommonTypes.o
endif

.PHONY: default lib debug gfortran aenet lammps clean help startUP modout finale removeObjects removeExec
