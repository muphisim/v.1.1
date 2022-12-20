include $(PETSC_DIR)/$(PETSC_ARCH)/lib/petsc/conf/petscvariables 
SHELL     := /bin/bash

################################################################
## Colordefinition
################################################################
NO_COLOR    = \x1b[0m
OK_COLOR    = \x1b[32;01m
WARN_COLOR  = \x1b[33;01m
ERROR_COLOR = \x1b[31;01m
CLEAN_COLOR = \x1b[94m 
TYPE :=SEQUENTIAL
METIS :=YES
CBLAS :=YES
MM :=NO
DEBUG := NO
# GENERAL FLAGS
CFLAGS= -Wall -std=c++0x  -frounding-math -Wl,-rpath,/usr/local/lib
FLAGS= -lgfortran -llapack -lblas -lm
PETSFLAGS=-lpetscsles -lpetscdm -lpetscmat -lpetscvec -lpetsc

ifeq ($(MM),YES)
MAXENT_LIB_PATH :=$(shell pwd)/lib/MAXENT-V1.4/lib
CFLAGS += -DMM -L$(MAXENT_LIB_PATH) -lmaxent
else
$(info The code is compiled without MM)
endif

ifeq ($(DEBUG),YES)
$(info The code is compiled in debug mode)	
CFLAGS += -g -D_GLIBCXX_DEBUG
else
$(info The code is compiled in optimum mode)
CFLAGS += -O3 -DNDEBUG -w
endif

# Folder paths
SOURCE_DIR=$(shell pwd)/src 
BUILD_DIR=$(shell pwd)/build
BIN_DIR=$(shell pwd)/bin
DOC_DIR=$(shell pwd)/doc
TEST_DIR=$(shell pwd)/applications

# Sequential and parallel versions
# Default values
ifeq ($(TYPE),SEQUENTIAL)
CC :=g++ # GNU Compiler
ifeq ($(PETSCMPI),YES)
CC=mpiCC # GNU Compiler
endif
else
ifeq ($(TYPE),PARALLEL)
CC=mpiCC # GNU Compiler
# Adding the flag for compilation time...
CFLAGS += -D$(TYPE)
ifeq ($(METIS),YES)
$(info The code is compiled with METIS)	
CFLAGS += -DMETIS -lmetis
else
$(info The code is not compiled with METIS, the mesh must be partitionned before running simulation)	
endif
else
$(error TYPE should be either SEQUENTIAL or PARALLEL. You sent $(TYPE))	
endif
endif

#FSI using preCICE
ifneq ($(PRECICE_DIR),)
$(info FSI is compiled using precice with PRECICE_DIR=$(PRECICE_DIR))	
CFLAGS += -DFSI -I$(PRECICE_DIR) -lprecice
endif

# Files type
vpath %.h $(SOURCE_DIR)
vpath %.cpp $(SOURCE_DIR)
vpath %.o $(BUILD_DIR)

# Add new constitutive models here to be compiled
CONS_OBJECTS= constitutiveModels.o stochasticHyperElasticNeoHookean.o stochasticHyperElasticStVenantKirchhoff.o hyperElasticStVenantKirchhoff.o hyperElasticStVenantKirchhoffThermoMechanics.o hyperElasticNeoHookean.o isoMorphogeneticGrowth.o fiberGrowth.o areaGrowth.o viscoElastic.o viscoElasticGrowth.o volMechDrivenGrowth.o viscoMechAxialGrowth.o viscoMechContractAxialGrowth.o FHN.o temperature.o 3DprintingThermoelasticityStVenantKirchhoff.o linearElastic.o constitutiveManager.o constitutiveAlgorithms.o J2plasticityFiniteStrain.o
CONS_SRC=$(CONS_OBJECTS:.o=.cpp)
CONS_DIR=$(shell pwd)/src/constitutiveModels

# Add new solver files here to be compiled
SOLV_OBJECTS= solvers.o solverExplicit.o nonlinearSystem.o solverImplicit.o boundaryConditions.o classNeumannBCs.o FSIAdapter.o dataExtraction.o timeStepping.o errorEstimation.o
SOLV_DIR=$(shell pwd)/src/solvers

# Add new spatial discretisation models here to be compiled
SPACE_OBJECTS= classSpatialPoint.o classNodes.o  GPsDistribution.o shapeFunctions.o classElements.o classGPs.o supportDomains.o classLinearTriangles.o classQuadTetrahedra.o classQuadTriangles.o classQuadLine.o classLinearTetrahedra.o classLinearLine.o classLinearSquare.o classLinearHexahedra.o classStates.o
SPACE_DIR=$(shell pwd)/src/spatialDiscretisation

# Add new input output files here to be compiled
IO_OBJECTS= commandLine.o inputFile.o output.o
IO_DIR=$(shell pwd)/src/IO


# Add new mathematical files here to be compiled
MATHS_OBJECTS= maths.o geometry.o
MATHS_DIR=$(shell pwd)/src/mathematics


# General objects
OBJECTS=configuration.o init.o 

OBJECTIVES=MuPhiSim
OBJECTS_MuPhiSim=MuPhiSim.o

# Name of the executable file
EXECUTABLE:=MuPhiSim

# General construction rules
# Default
all: print build_dirs $(OBJECTIVES)

print:
	@echo -e "$(OK_COLOR)MuPhiSim is being compiled in its $(WARN_COLOR)$(TYPE) $(OK_COLOR)version. TYPE=$(WARN_COLOR)SEQUENTIAL/PARALLEL$(OK_COLOR) when making"

build_dirs:
	@mkdir -p $(BUILD_DIR) $(BIN_DIR)
    
MuPhiSim: $(OBJECTS) $(SPACE_OBJECTS) $(IO_OBJECTS) $(CONS_OBJECTS) $(SOLV_OBJECTS) $(MATHS_OBJECTS) $(OBJECTS_MuPhiSim)
	@echo -e "$(OK_COLOR)Compiling and linking $(WARN_COLOR)$(EXECUTABLE) $(NO_COLOR)"
	@echo -e "$(OK_COLOR)The executable has been created and placed in the bin folder"
	@cd $(BUILD_DIR); $(CC) $^ $(CFLAGS) ${PETSC_LIB} $(LDFLAGS)  $(FLAGS) -o $(BIN_DIR)/$(EXECUTABLE)
# Prevous objects. Common
%.o: %.cpp %.h
	@echo -e "$(OK_COLOR)Compiling $< $(WARN_COLOR)$@$(NO_COLOR)"
	@$(CC) $(CFLAGS) $(FLAGS) -I$(SOURCE_DIR) -I$(CONS_DIR) -I$(MATHS_DIR) -I$(SOLV_DIR) -I$(SPACE_DIR) -I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(IO_DIR) -c $< -o $(BUILD_DIR)/$@
$(MATHS_OBJECTS):%.o: $(MATHS_DIR)/%.cpp $(MATHS_DIR)/%.h
	@echo -e "$(OK_COLOR)Compiling $< $(WARN_COLOR)$@$(NO_COLOR)"
	@$(CC) $(CFLAGS) $(FLAGS) -I$(SOURCE_DIR) -I$(CONS_DIR) -I$(MATHS_DIR) -I$(SOLV_DIR) -I$(SPACE_DIR) -I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(IO_DIR) -c $< -o $(BUILD_DIR)/$@
$(IO_OBJECTS):%.o: $(IO_DIR)/%.cpp $(IO_DIR)/%.h 
	@echo -e "$(OK_COLOR)Compiling $< $(WARN_COLOR)$@$(NO_COLOR)"
	@$(CC) $(CFLAGS)  $(FLAGS)  -I$(SOURCE_DIR) -I$(CONS_DIR) -I$(MATHS_DIR)  -I$(SOLV_DIR) -I$(SPACE_DIR) -I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(IO_DIR) -c $< -o $(BUILD_DIR)/$@
$(SPACE_OBJECTS):%.o: $(SPACE_DIR)/%.cpp $(SPACE_DIR)/%.h
	@echo -e "$(OK_COLOR)Compiling $< $(WARN_COLOR)$@$(NO_COLOR)"
	@$(CC) $(CFLAGS)  $(FLAGS) -I$(SOURCE_DIR) -I$(CONS_DIR) -I$(MATHS_DIR)  -I$(SOLV_DIR) -I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(SPACE_DIR) -I$(IO_DIR) -c $< -o $(BUILD_DIR)/$@
$(CONS_OBJECTS):%.o: $(CONS_DIR)/%.cpp $(CONS_DIR)/%.h
	@echo -e "$(OK_COLOR)Compiling $< $(WARN_COLOR)$@$(NO_COLOR)"
	@$(CC) $(CFLAGS)  $(FLAGS) -I$(SOURCE_DIR) -I$(CONS_DIR) -I$(MATHS_DIR)  -I$(SPACE_DIR) -I$(IO_DIR) -I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include  -c $< -o $(BUILD_DIR)/$@
$(SOLV_OBJECTS):%.o: $(SOLV_DIR)/%.cpp $(SOLV_DIR)/%.h
	@echo -e "$(OK_COLOR)Compiling $< $(WARN_COLOR)$@$(NO_COLOR)"
	@$(CC) $(CFLAGS)  $(FLAGS) -I$(SOURCE_DIR) -I$(CONS_DIR) -I$(SOLV_DIR) -I$(MATHS_DIR) -I$(SPACE_DIR) -I$(IO_DIR) -I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include -c $< -o $(BUILD_DIR)/$@
MuPhiSim.o: main.cpp
	@echo -e "$(OK_COLOR)Compiling $< $(WARN_COLOR)$@$(NO_COLOR)"
	@$(CC) $(CFLAGS)  $(FLAGS) -I$(SOURCE_DIR) -I$(SOLV_DIR) -I$(MATHS_DIR)  -I$(CONS_DIR) -I$(SPACE_DIR) -I$(IO_DIR) -I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include -c $< -o $(BUILD_DIR)/$@

# Cleaning rules
clean-all: clean-out clean-doc clean
# Cleaning rules for the documentation
clean-doc: 
	@echo -e "$(CLEAN_COLOR)"
	-cd $(DOC_DIR); rm -rf html latex
	@echo -e "$(NO_COLOR)" 
# General cleaning rules
clean:
	@echo -e "$(CLEAN_COLOR)"	
	-rm -rf $(BUILD_DIR)
	-rm -rf $(BIN_DIR)
	-cd $(SOURCE_DIR); rm -f *~
	-cd $(DOC_DIR); rm -rf html latex
	-rm -rf $(TEST_DIR)/input
	-rm -rf $(TEST_DIR)/output
	-cd $(IO_DIR); rm -f *~
	-cd $(SOLV_DIR); rm -f *~
	-cd $(SPACE_DIR); rm -f *~
	-cd $(CONS_DIR); rm -f *~
	-cd $(MATHS_DIR); rm -f *~
	@echo -e "$(NO_COLOR)"
# clean others but keep bin folder
clean-build:
	@echo -e "$(CLEAN_COLOR)"	
	-rm -rf $(BUILD_DIR)
	-cd $(SOURCE_DIR); rm -f *~
	-cd $(DOC_DIR); rm -rf html latex
	-rm -rf $(TEST_DIR)/input
	-rm -rf $(TEST_DIR)/output
	-cd $(IO_DIR); rm -f *~
	-cd $(SOLV_DIR); rm -f *~
	-cd $(SPACE_DIR); rm -f *~
	-cd $(CONS_DIR); rm -f *~
	-cd $(MATHS_DIR); rm -f *~
	@echo -e "$(NO_COLOR)"
# Cleaning rules for the documentation
clean-out:
	@echo -e "$(CLEAN_COLOR)" 
	-rm -rf $(TEST_DIR)/output
	@echo -e "$(NO_COLOR)" 
# Documentation rules
documentation:
	@echo -e "$(CLEAN_COLOR)"
	@cd $(DOC_DIR); doxygen
	@echo -e "$(OK_COLOR)The documentation has been created and placed in the doc folder"
	@echo -e "$(NO_COLOR)" 
ifeq ($(TYPE),SEQUENTIAL)
test:
	@echo -e "$(OK_COLOR)randomly take 2 tests and run in sequential"
	@mkdir -p ${TEST_DIR}/input ${TEST_DIR}/output 
	@for filename in $(shell cd ${TEST_DIR}/examples; find example*.inp -type f | shuf -n 2); do \
	    echo $${filename}; \
	    echo -e "$(OK_COLOR)Run test: $ $${filename}"; \
	    cd $(TEST_DIR); cp examples/$${filename} $(TEST_DIR)/input; $(BIN_DIR)/$(EXECUTABLE) -input $${filename};\
	done
else 
NCPUS :=4
test:
	@echo -e "$(OK_COLOR)randomly take 2 tests and run in parallel"
	@mkdir -p ${TEST_DIR}/input ${TEST_DIR}/output
	@for filename in $(shell cd ${TEST_DIR}/examples; find example*.inp -type f | shuf -n 2); do \
	    echo $${filename}; \
	    echo -e "$(OK_COLOR)Run test: $ $${filename}"; \
	    cd $(TEST_DIR); cp examples/$${filename} $(TEST_DIR)/input; mpirun -np $(NCPUS) $(BIN_DIR)/$(EXECUTABLE) -input $${filename};\
	done
endif


