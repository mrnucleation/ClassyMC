#SHELL = /bin/sh
CUR_DIR := $(shell pwd)
# ====================================
#        Compiler Options
# ====================================
#FC := mpif90
FC := /opt/openmpi/bin/mpif90
#FC := mpifort
#FC := gfortran
CC := mpicc
OPTIMIZE_FLAGS := -O3
#OPTIMIZE_FLAGS += -xHost
#OPTIMIZE_FLAGS += -ipo
#OPTIMIZE_FLAGS += -no-prec-div
#OPTIMIZE_FLAGS += -prof-gen -prof-dir=$(CUR_DIR)/profiling
#OPTIMIZE_FLAGS += -prof-use -prof-dir=$(CUR_DIR)/profiling
#DETAILEDDEBUG:= -fbacktrace -fcheck=all -g -ffree-line-length-0 -Og
DETAILEDDEBUG:= -check all -traceback -g -fpe3 -Og
#DEBUGFLAGS:= -fbacktrace -fcheck=all -g
#DEBUGFLAGS += -heap-arrays 1024
#DEBUGFLAGS += $(DETAILEDDEBUG)
#DEBUGFLAGS += -fpe3
#DEBUGFLAGS += -pg 
#DEBUGFLAGS += -ffpe-trap=invalid
#DEBUGFLAGS += -Wunused-parameter 
#DEBUGFLAGS := -fimplicit-none  -Wline-truncation -Wcharacter-truncation -Wsurprising -Waliasing -fwhole-file -fcheck=all -fbacktrace
COMPFLAGS := $(DEBUGFLAGS) $(OPTIMIZE_FLAGS)


# ====================================
#        Directory List
# ====================================

SRC := $(CUR_DIR)/src
OBJ := $(CUR_DIR)/objects


# Define file extensions
.SUFFIXES:
.SUFFIXES: .f .f90 .o .mod 
# ====================================
#        Compiler specific commands
# ====================================

#MODFLAGS := -I $(MODS) -J $(MODS)
# ====================================
#        Source Files
# ====================================

SRC_MAIN := $(SRC)/Common.f90\
        		$(SRC)/Common_BoxData.f90\
        		$(SRC)/Common_ECalc.f90\
        		$(SRC)/Common_MolDef.f90\
        		$(SRC)/Common_Sampling.f90\
        		$(SRC)/ConstraintClass.f90\
        		$(SRC)/DistanceCriteria.f90\
        		$(SRC)/AcceptRule.f90\
        		$(SRC)/Metropolis.f90\
        		$(SRC)/Main.f90\
        		$(SRC)/SimpleBox.f90\
        		$(SRC)/CubicBox.f90\
        		$(SRC)/RandomNew.f90\
        		$(SRC)/FF_LJ_Cut.f90\
        		$(SRC)/MoveClass.f90\
        		$(SRC)/AtomTranslation.f90\
        		$(SRC)/VariablePrecision.f90\
 	        	$(SRC)/Script_Forcefield.f90\
	        	$(SRC)/Script_Main.f90\
	        	$(SRC)/Script_SimBoxes.f90\
            	$(SRC)/Input_Format.f90\
 	        	$(SRC)/NeighList.f90

SRC_TEMPLATE := $(SRC)/Template_SimBox.f90\
        		$(SRC)/Template_Forcefield.f90

SRC_COMPLETE := $(SRC_TEMPLATE) $(SRC_MAIN) 

# ====================================
#        Object Files
# ====================================
OBJ_MAIN:=$(patsubst $(SRC)/%.f90, $(OBJ)/%.o, $(SRC_MAIN))

OBJ_TEMPLATE:=$(patsubst $(SRC)/%.f90, $(OBJ)/%.o, $(SRC_TEMPLATE))

OBJ_COMPLETE:= ${OBJ_TEMPLATE} $(OBJ_MAIN) 
# ====================================
#        Compile Commands
# ====================================


.f90.o :     
		@echo Creating $<
		@$(FC) $(COMPFLAGS) $(MODFLAGS) -c -o $@ $<


.f90.mod :     
		@echo Creating $<
		@$(FC) $(COMPFLAGS) $(MODFLAGS) -c -o $@ $<

$(OBJ)/%.o: $(SRC)/%.f90
		@echo Creating $<
		@$(FC) $(COMPFLAGS) $(MODFLAGS) -c -o $@ $<

$(OBJ)/%.o: $(TEMPLATE)/%.f90
		@echo Creating $<
		@$(FC) $(COMPFLAGS) $(MODFLAGS) -c -o $@ $<

# ====================================
#        Compile Commands
# ====================================
default: startUP classyMC finale
debug: startUP_debug classyMC_debug finale
neat: startUP classyMC removeObject finale
clean: removeObjects removeExec finale    
       
classyMC: $(OBJ_COMPLETE) 
		@echo =============================================
		@echo     Compiling and Linking Source Files
		@echo =============================================	
		@$(FC) $(COMPFLAGS) $(MODFLAGS)  $^ -o $@ 	
	

classyMC_debug: $(OBJ_COMPLETE) 
		@echo =============================================
		@echo     Compiling and Linking Source Files
		@echo =============================================	
		@$(FC) $(DETAILEDDEBUG) $(MODFLAGS)  $^ -o $@ 	
		
startUP:
		@echo ==================================================================
		@echo ---------------------- Begin ---------------------------------		
		@echo Current Directory:$(CUR_DIR)		
		@echo Compiler and Flags used:	$(FC) $(COMPFLAGS) 		
		@echo		

startUP_debug:
		@echo ==================================================================
		@echo ---------------------- Begin ---------------------------------		
		@echo Current Directory:$(CUR_DIR)		
		@echo Compiler and Flags used:	$(FC) $(DETAILEDDEBUG) 		
		@echo		

finale:
		@echo
		@echo ---------------------- Finished! ---------------------------------
		@echo ==================================================================		
     
removeObjects:
		@echo =============================================
		@echo            Cleaning Directory
		@echo =============================================		
		@echo		
		@rm -f ./*.o ./*.mod				
		@rm -f $(SRC)/*.o $(SRC)/*.mod		
		@rm -f $(MODS)/*.o $(MODS)/*.mod			
		@rm -f $(SRC)/*/*.o $(SRC)/*/*.mod
		@rm -f $(OBJ)/*.o		

removeExec:
		@rm -f $(CUR_DIR)/classyMC
		@rm -f $(CUR_DIR)/classyMC_debug
		@rm -f $(CUR_DIR)/classyMC.exe


# ====================================
#        Dependencies
# ====================================
$(OBJ)/Common.o: $(OBJ)/VariablePrecision.o
$(OBJ)/Common_BoxData.o: $(OBJ)/SimpleBox.o $(OBJ)/ConstraintClass.o
$(OBJ)/Common_ECalc.o: $(OBJ)/Template_Forcefield.o $(OBJ)/Common.o

$(OBJ)/Template_SimBox.o: $(OBJ)/Common.o ${OBJ}/Input_Format.o

$(OBJ)/SimpleBox.o: $(OBJ)/Common.o $(OBJ)/NeighList.o $(OBJ)/Input_Format.o $(OBJ)/Common_ECalc.o $(OBJ)/Template_SimBox.o $(OBJ)/ConstraintClass.o
$(OBJ)/CubicBox.o: $(OBJ)/SimpleBox.o

$(OBJ)/ConstraintClass.o: ${OBJ}/Template_SimBox.o 

$(OBJ)/Main.o: $(OBJ)/Common.o  $(OBJ)/Units.o  $(OBJ)/Script_Main.o $(OBJ)/AtomTranslation.o $(OBJ)/RandomNew.o
$(OBJ)/Template_Forcefield.o: $(OBJ)/Common.o  $(OBJ)/Common_MolDef.o $(OBJ)/Template_SimBox.o

$(OBJ)/AtomTranslation.o: $(OBJ)/Common.o $(OBJ)/Common_BoxData.o $(OBJ)/SimpleBox.o $(OBJ)/RandomNew.o $(OBJ)/MoveClass.o $(OBJ)/ConstraintClass.o

$(OBJ)/Script_Main.o: $(OBJ)/Common_BoxData.o $(OBJ)/Script_Forcefield.o $(OBJ)/CubicBox.o $(OBJ)/Script_SimBoxes.o
$(OBJ)/Script_Forcefield.o: ${OBJ}/Input_Format.o ${OBJ}/Template_Forcefield.o ${OBJ}/FF_LJ_Cut.o

$(OBJ)/RandomNew.o: $(OBJ)/Common.o
$(OBJ)/Metropolis.o: $(OBJ)/RandomNew.o
$(OBJ)/Common_Sampling.o: $(OBJ)/AcceptRule.o $(OBJ)/Metropolis.o
