#SHELL = /bin/sh
CUR_DIR := $(shell pwd)
# ====================================
#        Compiler Options
# ====================================
FC := mpif90
#FC := /opt/openmpi/bin/mpif90
#FC := mpifort
#FC := gfortran
CC := mpicc
OPTIMIZE_FLAGS := -O3
#OPTIMIZE_FLAGS += -xHost
#OPTIMIZE_FLAGS += -ipo
#OPTIMIZE_FLAGS += -no-prec-div
#OPTIMIZE_FLAGS += -prof-gen -prof-dir=$(CUR_DIR)/profiling
#OPTIMIZE_FLAGS += -prof-use -prof-dir=$(CUR_DIR)/profiling
DETAILEDDEBUG:= -fbacktrace -fcheck=all -g -ffree-line-length-0 -Og
#DETAILEDDEBUG:= -check all -traceback -g -fpe3 -Og
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
        		$(SRC)/Common_TrajData.f90\
        		$(SRC)/Common_ECalc.f90\
        		$(SRC)/Common_MolDef.f90\
        		$(SRC)/Common_Sampling.f90\
        		$(SRC)/Common_MCMoves.f90\
        		$(SRC)/Common_NeighList.f90\
         		$(SRC)/Constrain_DistCriteria.f90\
        		$(SRC)/Sampling_Metropolis.f90\
        		$(SRC)/Move_AtomTranslation.f90\
        		$(SRC)/ExeptionHandling.f90\
        		$(SRC)/Box_SimpleBox.f90\
        		$(SRC)/Box_CubicBox.f90\
        		$(SRC)/Box_OrthoBox.f90\
        		$(SRC)/Box_Ultility.f90\
        		$(SRC)/RandomNew.f90\
        		$(SRC)/Constrain_HardWall.f90\
        		$(SRC)/FF_LJ_Cut.f90\
        		$(SRC)/FF_LJ_Cut_NoNei.f90\
        		$(SRC)/FF_Tersoff.f90\
 	        	$(SRC)/Script_Forcefield.f90\
 	        	$(SRC)/Script_FieldType.f90\
 	        	$(SRC)/Script_NeighType.f90\
 	        	$(SRC)/Script_TrajType.f90\
 	        	$(SRC)/Script_LoadCoords.f90\
	        	$(SRC)/Script_Main.f90\
	        	$(SRC)/Script_Sampling.f90\
	        	$(SRC)/Script_Initialize.f90\
	        	$(SRC)/Script_MCMoves.f90\
	        	$(SRC)/Script_SimBoxes.f90\
	        	$(SRC)/Output_DumpCoords.f90\
	        	$(SRC)/Traj_XYZFormat.f90\
						$(SRC)/Input_Format.f90\
 	        	$(SRC)/Neigh_RSqList.f90\
        		$(SRC)/VariablePrecision.f90\
        		$(SRC)/Main.f90\
        		$(SRC)/Units.f90

SRC_TEMPLATE := $(SRC)/Template_SimBox.f90\
                $(SRC)/Template_Constraint.f90\
	              $(SRC)/Template_Forcefield.f90\
								$(SRC)/Template_AcceptRule.f90\
								$(SRC)/Template_NeighList.f90\
								$(SRC)/Template_Trajectory.f90\
                $(SRC)/Template_MultiBoxMove.f90\
								$(SRC)/Template_MoveClass.f90

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
$(OBJ)/Common_BoxData.o: $(OBJ)/Box_SimpleBox.o 
$(OBJ)/Common_ECalc.o: $(OBJ)/Template_Forcefield.o $(OBJ)/Common.o
$(OBJ)/Common_Sampling.o: $(OBJ)/Template_AcceptRule.o $(OBJ)/Sampling_Metropolis.o

$(OBJ)/Neigh_RSqList.o: $(OBJ)/Common_BoxData.o $(OBJ)/Template_NeighList.o $(OBJ)/Common_NeighList.o

$(OBJ)/Template_Constraint.o: $(OBJ)/Template_SimBox.o 
$(OBJ)/Template_SimBox.o: $(OBJ)/Common.o ${OBJ}/Input_Format.o $(OBJ)/Template_NeighList.o
$(OBJ)/Template_MultiBoxMove.o: $(OBJ)/Template_MoveClass.o
$(OBJ)/Template_MoveClass.o: $(OBJ)/Common.o ${OBJ}/Box_SimpleBox.o
$(OBJ)/Template_Forcefield.o: $(OBJ)/Common.o  $(OBJ)/Common_MolDef.o $(OBJ)/Template_SimBox.o

$(OBJ)/Box_SimpleBox.o: $(OBJ)/Common.o $(OBJ)/Template_NeighList.o $(OBJ)/Input_Format.o $(OBJ)/Common_ECalc.o $(OBJ)/Template_SimBox.o $(OBJ)/Template_Constraint.o 
$(OBJ)/Box_CubicBox.o: $(OBJ)/Box_SimpleBox.o
$(OBJ)/Box_OrthoBox.o: $(OBJ)/Box_SimpleBox.o
$(OBJ)/Box_Utility.o: $(OBJ)/Box_SimpleBox.o


$(OBJ)/Move_AtomTranslation.o: $(OBJ)/Common.o $(OBJ)/Common_BoxData.o $(OBJ)/Box_SimpleBox.o $(OBJ)/RandomNew.o $(OBJ)/Template_MoveClass.o $(OBJ)/Template_Constraint.o $(OBJ)/Box_Ultility.o

$(OBJ)/Script_Main.o: $(OBJ)/Units.o $(OBJ)/Common_BoxData.o $(OBJ)/Script_Forcefield.o $(OBJ)/Box_CubicBox.o $(OBJ)/Script_SimBoxes.o $(OBJ)/Script_Sampling.o $(OBJ)/Script_MCMoves.o $(OBJ)/Script_Initialize.o
$(OBJ)/Script_Forcefield.o: ${OBJ}/Input_Format.o ${OBJ}/Template_Forcefield.o  ${OBJ}/Move_AtomTranslation.o ${OBJ}/Units.o $(OBJ)/Script_FieldType.o
$(OBJ)/Script_LoadCoords.o: ${OBJ}/Script_SimBoxes.o
$(OBJ)/Script_FieldType.o: ${OBJ}/Input_Format.o ${OBJ}/Template_Forcefield.o ${OBJ}/FF_LJ_Cut.o ${OBJ}/FF_LJ_Cut_NoNei.o ${OBJ}/FF_Tersoff.o ${OBJ}/Move_AtomTranslation.o $(OBJ)/Common_ECalc.o
$(OBJ)/Script_TrajType.o: ${OBJ}/Common_TrajData.o ${OBJ}/Template_Trajectory.o ${OBJ}/Traj_XYZFormat.o 
$(OBJ)/Script_NeighType.o: ${OBJ}/Neigh_RSqList.o $(OBJ)/Common_BoxData.o

$(OBJ)/RandomNew.o: $(OBJ)/Common.o

$(OBJ)/Sampling_Metropolis.o: $(OBJ)/RandomNew.o

$(OBJ)/Main.o: $(OBJ)/Common.o  $(OBJ)/Units.o  $(OBJ)/Script_Main.o $(OBJ)/Move_AtomTranslation.o $(OBJ)/RandomNew.o 


