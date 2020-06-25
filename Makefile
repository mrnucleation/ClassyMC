#SHELL = /bin/sh
CUR_DIR := $(shell pwd)
# ====================================
#        Compiler Options
# ====================================
FC := mpif90
#FC := /opt/openmpi/bin/mpif90
#FC := mpifort
#FC := gfortran
AR := ar
CC := mpicc
OPTIMIZE_FLAGS_IFORT := -O3
OPTIMIZE_FLAGS_IFORT += -xHost
#OPTIMIZE_FLAGS_IFORT += -ipo
OPTIMIZE_FLAGS_IFORT += -no-prec-div
OPTIMIZE_FLAGS_IFORT += -no-wrap-margin
OPTIMIZE_FLAGS_IFORT += -fpp
#OPTIMIZE_FLAGS_IFORT += -fpe0
#OPTIMIZE_FLAGS_IFORT += -pg
OPTIMIZE_FLAGS_IFORT += -traceback
#OPTIMIZE_FLAGS_IFORT += -prof-gen -prof-dir=$(CUR_DIR)/profiling
#OPTIMIZE_FLAGS_IFORT += -prof-use -prof-dir=$(CUR_DIR)/profiling

OPTIMIZE_FLAGS_GFORT := -O3 -cpp -g
OPTIMIZE_FLAGS_GFORT += -fbacktrace -fcheck=bounds -ffree-line-length-512
OPTIMIZE_FLAGS_GFORT += -ffpe-trap=overflow,invalid,zero
#OPTIMIZE_FLAGS_GFORT += -pg
#OPTIMIZE_FLAGS_GFORT += -lblas -llapack

LIBRARY_FLAGS := -shared -fpic

DETAILEDDEBUG_GFORT:= -fbacktrace -fcheck=all -g -ffree-line-length-0 -Og -cpp -ffpe-trap=overflow,invalid,zero 
DETAILEDDEBUG_IFORT:= -check all -traceback -g -fpe0 -O0 -fp-stack-check -debug all -ftrapuv -fpp -no-wrap-margin
#DEBUGFLAGS:= -check all -warn -traceback -g -fpe0 -O0 -fp-stack-check -debug all -ftrapuv 
#DEBUGFLAGS:= -fbacktrace -fcheck=all -g
#DEBUGFLAGS += -fpe0
#DEBUGFLAGS += -pg 
#DEBUGFLAGS += -ffpe-trap=invalid
#DEBUGFLAGS += -Wunused-parameter 
#DEBUGFLAGS := -fimplicit-none  -Wline-truncation -Wcharacter-truncation -Wsurprising -Waliasing -fwhole-file -fcheck=all -fbacktrace
#COMPFLAGS := $(DEBUGFLAGS) $(OPTIMIZE_FLAGS)


PACKAGE_FLAGS := -DPARALLEL

# ====================================
#        Directory List
# ====================================

SRC := $(CUR_DIR)/src
LIB := $(CUR_DIR)/lib
MOD := $(CUR_DIR)/mods
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
        		$(SRC)/C_To_Fotran.f90\
        		$(SRC)/Common_BoxData.f90\
        		$(SRC)/Common_TrajData.f90\
        		$(SRC)/Common_Analysis.f90\
        		$(SRC)/Common_ECalc.f90\
        		$(SRC)/Common_MolDef.f90\
        		$(SRC)/Common_Sampling.f90\
        		$(SRC)/Common_MCMoves.f90\
        		$(SRC)/Common_NeighList.f90\
        		$(SRC)/Debug.f90\
         		$(SRC)/Constrain_DistCriteria.f90\
         		$(SRC)/Constrain_EnergyCeiling.f90\
         		$(SRC)/Constrain_EnergyFloor.f90\
         		$(SRC)/Constrain_FreezeType.f90\
         		$(SRC)/Constrain_HardWall.f90\
	        	$(SRC)/SearchSort.f90\
        		$(SRC)/Sampling_AcceptAll.f90\
        		$(SRC)/Sampling_AcceptNone.f90\
        		$(SRC)/Sampling_Metropolis.f90\
        		$(SRC)/Sampling_MinMetrop.f90\
        		$(SRC)/Sampling_Nested.f90\
        		$(SRC)/Sampling_Umbrella.f90\
        		$(SRC)/Sampling_UmbrellaWHAM.f90\
        		$(SRC)/Move_MC_AVBMC.f90\
        		$(SRC)/Move_MC_AnisoVol.f90\
        		$(SRC)/Move_MC_AtomExchange.f90\
        		$(SRC)/Move_MC_AtomTranslation.f90\
        		$(SRC)/Move_MC_ParticleExchange.f90\
        		$(SRC)/Move_MC_BasicSwap.f90\
        		$(SRC)/Move_MC_Delete.f90\
        		$(SRC)/Move_MC_MolTranslation.f90\
        		$(SRC)/Move_MC_IsoVol.f90\
        		$(SRC)/Move_MC_PlaneRotate.f90\
        		$(SRC)/Move_MC_PlaneTranslate.f90\
        		$(SRC)/Move_MC_ThermoLambda.f90\
        		$(SRC)/Move_MC_UBSwap.f90\
        		$(SRC)/Move_MC_VolExchange.f90\
        		$(SRC)/ExeptionHandling.f90\
        		$(SRC)/MolSearch.f90\
        		$(SRC)/Analysis_BlockAverage.f90\
        		$(SRC)/Analysis_DensityOfStates.f90\
        		$(SRC)/Analysis_ClusterSize.f90\
        		$(SRC)/Analysis_DistPair.f90\
        		$(SRC)/Analysis_RDF.f90\
        		$(SRC)/Analysis_ThermoAverage.f90\
        		$(SRC)/Analysis_ThermoIntegration.f90\
        		$(SRC)/Box_Presets.f90\
        		$(SRC)/Box_SimpleBox.f90\
        		$(SRC)/Box_CubicBox.f90\
        		$(SRC)/Box_OrthoBox.f90\
        		$(SRC)/Box_Ultility.f90\
        		$(SRC)/RandomNew.f90\
        		$(SRC)/FF_AENet.f90\
        		$(SRC)/FF_Einstein.f90\
        		$(SRC)/FF_HardSphere.f90\
        		$(SRC)/FF_Hybrid.f90\
        		$(SRC)/FF_LJ_Cut.f90\
        		$(SRC)/FF_LJWall.f90\
        		$(SRC)/FF_LJ_Shift.f90\
        		$(SRC)/FF_LJ_Ele_Cut.f90\
        		$(SRC)/FF_Pedone.f90\
        		$(SRC)/FF_Tersoff.f90\
        		$(SRC)/FF_ThermoInt.f90\
        		$(SRC)/Intra_AngleRidgid.f90\
        		$(SRC)/Intra_BondRidgid.f90\
        		$(SRC)/Intra_BondHarmonic.f90\
        		$(SRC)/Intra_TorsionRidgid.f90\
        		$(SRC)/MolCon_LinearCBMC.f90\
        		$(SRC)/MolCon_RidgidRegrowth.f90\
        		$(SRC)/MolCon_SimpleRegrowth.f90\
 	        	$(SRC)/Script_AnalysisType.f90\
 	        	$(SRC)/Script_AngleType.f90\
 	        	$(SRC)/Script_BondType.f90\
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
 	        	$(SRC)/Neigh_RSqList.f90\
        		$(SRC)/VariablePrecision.f90\
        		$(SRC)/Sim_Minimize.f90\
        		$(SRC)/Sim_MonteCarlo.f90\
        		$(SRC)/Sim_GeneticAlgor.f90\
        		$(SRC)/Sim_LibControl.f90\
        		$(SRC)/Units.f90

SRC_TEMPLATE := $(SRC)/Template_Master.f90\
				$(SRC)/Template_AcceptRule.f90\
				$(SRC)/Template_Analysis.f90\
				$(SRC)/Template_BondFF.f90\
				$(SRC)/Template_AngleFF.f90\
				$(SRC)/Template_TorsionFF.f90\
				$(SRC)/Template_Constraint.f90\
				$(SRC)/Template_Forcefield.f90\
	            $(SRC)/Template_SimBox.f90\
				$(SRC)/Template_IntraFF.f90\
				$(SRC)/Template_NeighList.f90\
				$(SRC)/Template_Trajectory.f90\
				$(SRC)/Template_MultiBoxMove.f90\
				$(SRC)/Template_MoveClass.f90\
				$(SRC)/Template_MolConstructor.f90


#SRC_COMPLETE := $(SRC_TEMPLATE) $(SRC_MAIN) 
SRC_COMPLETE := $(SRC_TEMPLATE) $(SRC_MAIN) 

# ====================================
#        Object Files
# ====================================
OBJ_MAIN:=$(patsubst $(SRC)/%.f90, $(OBJ)/%.o, $(SRC_MAIN))
OBJ_TEMPLATE:=$(patsubst $(SRC)/%.f90, $(OBJ)/%.o, $(SRC_TEMPLATE))
OBJ_AENET:=$(patsubst $(SRC)/%.f90, $(OBJ)/%.o, $(SRC_AENET))
OBJ_COMPLETE:= $(OBJ_TEMPLATE) $(OBJ_MAIN) 
# ====================================
#        Compile Commands
# ====================================
default: COMPFLAGS := $(OPTIMIZE_FLAGS_IFORT) $(PACKAGE_FLAGS)
default: COMPFLAGS += $(DEBUGFLAGS)
default: SRC_COMPLETE += $(SRC)/Main.f90
default: startUP  classyMC modout finale 

aenet: COMPFLAGS := $(OPTIMIZE_FLAGS_IFORT) $(PACKAGE_FLAGS)
#aenet: COMPFLAGS := $(OPTIMIZE_FLAGS_GFORT) $(PACKAGE_FLAGS)
aenet: COMPFLAGS += $(DEBUGFLAGS)
aenet: COMPFLAGS += -DAENET -static-libgfortran -llapack -lblas
aenet: startUP  classyMCAENet  modout finale 

#lib: COMPFLAGS := $(OPTIMIZE_FLAGS_GFORT) $(LIBRARY_FLAGS)
lib: COMPFLAGS := $(OPTIMIZE_FLAGS_IFORT) $(LIBRARY_FLAGS)
lib:  startUP  libclassymc.so  modout finale 

lib_debug: COMPFLAGS := $(DETAILEDDEBUG_IFORT) $(LIBRARY_FLAGS)
lib_debug: startUP  libclassymc.so  modout finale

gfortran: COMPFLAGS := $(OPTIMIZE_FLAGS_GFORT) $(PACKAGE_FLAGS)
gfortran: COMPFLAGS += $(DEBUGFLAGS)
gfortran:  startUP classyMC  modout finale 

debug: COMPFLAGS := $(DETAILEDDEBUG_IFORT) $(PACKAGE_FLAGS)
debug: startUP_debug classyMC_debug  modout finale

#debugaenet: COMPFLAGS := $(DETAILEDDEBUG_GFORT)
debugaenet: COMPFLAGS := $(DETAILEDDEBUG_IFORT) $(PACKAGE_FLAGS)
debugaenet: COMPFLAGS += -DAENET
debugaenet: startUP classyMCAENet  modout finale


debug_gfortran: COMPFLAGS := $(DETAILEDDEBUG_GFORT) $(PACKAGE_FLAGS)
debug_gfortran: startUP_debug classyMC_debug  modout finale

#neat: startUP classyMC removeObject finale

clean:  removeObjects removeExec finale    

# ====================================
#        Compile Commands
# ====================================

%.o: %.mod

.f90.o :     
		@echo Creating $<
		@$(FC) $(COMPFLAGS) $(MODFLAGS) -c -o $@ $<


.f90.mod :     
		@echo Creating $<
		@$(FC) $(COMPFLAGS) $(MODFLAGS) -c -o $@ $<

$(OBJ)/%.o: $(SRC)/%.f90
		@echo Creating $<
		@$(FC) $(COMPFLAGS) $(MODFLAGS) -c $< -o $@ 

$(OBJ)/%.o: $(TEMPLATE)/%.f90
		@echo Creating $<
		@$(FC) $(COMPFLAGS) $(MODFLAGS) -c $< -o $@ 


libclassymc.so: $(OBJ_COMPLETE) 
		@echo =============================================
		@echo     Compiling and Linking Source Files
		@echo =============================================	
		@$(FC) $(COMPFLAGS) $(MODFLAGS)  $^ -o $@ 	

       
classyMC: $(OBJ_COMPLETE) $(SRC)/Main.f90 
		@echo =============================================
		@echo     Compiling and Linking Source Files
		@echo =============================================	
		@$(FC) $(COMPFLAGS) $(MODFLAGS)  $^ -o $@ 	
	
classyMCAENet: $(OBJ_COMPLETE) $(SRC)/Main.f90 $(LIB)/libaenet.a
		@echo =============================================
		@echo     Compiling and Linking Source Files
		@echo =============================================	
		@$(FC) $(COMPFLAGS) $(MODFLAGS)  $^ -o $@ 	
	
classyMC_debug: $(OBJ_COMPLETE) $(SRC)/Main.f90
	    
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
		@mv $(MOD)/*.mod $(CUR_DIR)/ || echo 

startUP_debug:
		@echo ==================================================================
		@echo ---------------------- Begin ---------------------------------		
		@echo Current Directory:$(CUR_DIR)		
		@echo Compiler and Flags used:	$(FC) $(COMPFLAGS)
		@mv $(MOD)/*.mod $(CUR_DIR)/ || echo 
		@echo		

modout:
		@mv  $(CUR_DIR)/*.mod $(MOD)/ || echo

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
		@rm -f $(MOD)/*.o $(MOD)/*.mod			
		@rm -f $(SRC)/*/*.o $(SRC)/*/*.mod
		@rm -f $(OBJ)/*.o		

removeExec:
		@rm -f $(CUR_DIR)/libclassymc.so
		@rm -f $(CUR_DIR)/classyMC
		@rm -f $(CUR_DIR)/classyMCAENet
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
$(OBJ)/Common_MolDef.o: $(OBJ)/Template_MolConstructor.o  $(OBJ)/Template_BondFF.o

$(OBJ)/Neigh_RSqList.o: $(OBJ)/Common_BoxData.o $(OBJ)/Template_NeighList.o $(OBJ)/Common_NeighList.o

$(OBJ)/Template_Constraint.o: $(OBJ)/Template_SimBox.o 
$(OBJ)/Template_AcceptRule.o: $(OBJ)/Common.o $(OBJ)/Input_Format.o $(OBJ)/Template_SimBox.o 
$(OBJ)/Template_Anaylsis.o: $(OBJ)/Box_SimpleBox.o
$(OBJ)/Template_SimBox.o: $(OBJ)/Common.o ${OBJ}/Input_Format.o $(OBJ)/Template_NeighList.o
$(OBJ)/Template_Master.o: $(OBJ)/VariablePrecision.o
$(OBJ)/Template_MultiBoxMove.o: $(OBJ)/Template_MoveClass.o
$(OBJ)/Template_MoveClass.o: $(OBJ)/Common.o ${OBJ}/Box_SimpleBox.o
$(OBJ)/Template_Forcefield.o: $(OBJ)/Common.o  $(OBJ)/Common_MolDef.o $(OBJ)/Template_SimBox.o
$(OBJ)/Template_NeighList.o: $(OBJ)/SearchSort.o
$(OBJ)/Template_MolConstructor.o: $(OBJ)/Template_SimBox.o 
$(OBJ)/Template_AngleFF.o: $(OBJ)/Template_IntraFF.o
$(OBJ)/Template_BondFF.o: $(OBJ)/Template_IntraFF.o
$(OBJ)/Template_TorsionFF.o: $(OBJ)/Template_IntraFF.o

$(OBJ)/Analysis_ThermoIntegration.o: $(OBJ)/FF_ThermoInt.o
$(OBJ)/Box_SimpleBox.o: $(OBJ)/Common.o $(OBJ)/Template_NeighList.o $(OBJ)/Input_Format.o $(OBJ)/Common_ECalc.o $(OBJ)/Template_SimBox.o $(OBJ)/Template_Constraint.o $(OBJ)/Units.o $(OBJ)/Common_NeighList.o
$(OBJ)/Box_CubicBox.o: $(OBJ)/Box_SimpleBox.o
$(OBJ)/Box_OrthoBox.o: $(OBJ)/Box_SimpleBox.o
$(OBJ)/Box_Utility.o: $(OBJ)/Box_SimpleBox.o
$(OBJ)/Box_Presets.o: $(OBJ)/Box_OrthoBox.o $(OBJ)/Box_CubicBox.o


$(OBJ)/Move_MC_AVBMC.o: $(OBJ)/Common.o $(OBJ)/Box_Ultility.o
$(OBJ)/Move_MC_AtomTranslation.o: $(OBJ)/Common.o $(OBJ)/Common_BoxData.o $(OBJ)/Box_SimpleBox.o $(OBJ)/RandomNew.o $(OBJ)/Template_MoveClass.o $(OBJ)/Template_Constraint.o $(OBJ)/Box_Ultility.o
$(OBJ)/Move_MC_IsoVol.o: $(OBJ)/Common.o $(OBJ)/Common_BoxData.o $(OBJ)/Box_CubicBox.o $(OBJ)/Box_OrthoBox.o $(OBJ)/RandomNew.o $(OBJ)/Template_MoveClass.o $(OBJ)/Template_Constraint.o $(OBJ)/Box_Ultility.o
$(OBJ)/Move_MC_AnisoVol.o: $(OBJ)/Common.o $(OBJ)/Common_BoxData.o $(OBJ)/Box_CubicBox.o $(OBJ)/Box_OrthoBox.o $(OBJ)/RandomNew.o $(OBJ)/Template_MoveClass.o $(OBJ)/Template_Constraint.o $(OBJ)/Box_Ultility.o
$(OBJ)/Move_MC_AtomExchange.o: $(OBJ)/Common.o $(OBJ)/Common_BoxData.o $(OBJ)/Box_SimpleBox.o $(OBJ)/RandomNew.o $(OBJ)/Template_MoveClass.o $(OBJ)/Box_Ultility.o
$(OBJ)/Move_MC_ThermoLambda.o: $(OBJ)/FF_ThermoInt.o $(OBJ)/Analysis_ThermoIntegration.o 
$(OBJ)/Move_GA_AtomExchange.o: $(OBJ)/Common.o $(OBJ)/Common_BoxData.o $(OBJ)/Box_Ultility.o

$(OBJ)/MolCon_SimpleRegrowth.o: $(OBJ)/Template_MolConstructor.o

$(OBJ)/Script_Main.o: $(OBJ)/Units.o $(OBJ)/Common_BoxData.o $(OBJ)/Script_Forcefield.o $(OBJ)/Box_CubicBox.o $(OBJ)/Script_SimBoxes.o $(OBJ)/Script_Sampling.o $(OBJ)/Script_MCMoves.o $(OBJ)/Script_Initialize.o $(OBJ)/Script_NeighType.o $(OBJ)/Script_TrajType.o $(OBJ)/Sim_MonteCarlo.o $(OBJ)/Sim_Minimize.o

$(OBJ)/Script_Forcefield.o: ${OBJ}/Input_Format.o ${OBJ}/Template_Forcefield.o  ${OBJ}/Move_MC_AtomTranslation.o ${OBJ}/Units.o $(OBJ)/Script_FieldType.o $(OBJ)/Script_BondType.o $(OBJ)/Script_AngleType.o $(OBJ)/Script_RegrowType.o 
$(OBJ)/Script_LoadCoords.o: ${OBJ}/Script_SimBoxes.o
$(OBJ)/Script_FieldType.o: ${OBJ}/Input_Format.o ${OBJ}/Template_Forcefield.o ${OBJ}/FF_LJ_Cut.o ${OBJ}/Move_MC_AtomTranslation.o $(OBJ)/Common_ECalc.o
$(OBJ)/Script_TrajType.o: ${OBJ}/Common_TrajData.o ${OBJ}/Template_Trajectory.o ${OBJ}/Traj_XSF.o ${OBJ}/Traj_XYZFormat.o $(OBJ)/Traj_LAMMPSDump.o $(OBJ)/Traj_POSCAR.o
$(OBJ)/Script_NeighType.o: ${OBJ}/Neigh_RSqList.o $(OBJ)/Neigh_CellRSqList.o $(OBJ)/Common_BoxData.o

$(OBJ)/RandomNew.o: $(OBJ)/Common.o $(OBJ)/Units.o

$(OBJ)/Sampling_Umbrella.o: $(OBJ)/Sampling_UmbrellaWHAM.o
$(OBJ)/Sampling_Metropolis.o: $(OBJ)/RandomNew.o

$(OBJ)/Main.o: $(OBJ)/Sim_MonteCarlo.o $(OBJ)/Sim_Minimize.o
$(OBJ)/Sim_Library.o: $(OBJ)/Script_Main.o


$(OBJ)/Sim_MonteCarlo.o: $(OBJ)/Common.o  $(OBJ)/Units.o  $(OBJ)/Move_MC_AtomTranslation.o $(OBJ)/RandomNew.o $(OBJ)/Common_TrajData.o $(OBJ)/Output_DumpCoords.o


