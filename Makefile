CXX = g++
CC  = gcc
FF  = gfortran
Opt = -Ofast

include Makefile.in

CXX_FLAGS = -ansi -fexceptions -DMATLAB_MEX_FILE -std=c++11 -fopenmp -march=native \
			-D_GNU_SOURCE -fPIC -fno-omit-frame-pointer -Wno-write-strings -pthread\
			$(Opt) -DNDEBUG -fopenmp -ffast-math
			
CXX_INCLUDE = -I./include/mexplus \
              -I./include/ \
              -I$(MATLAB_ROOT)extern/include \
              -I$(MATLAB_ROOT)simulink/include
              
MATLAB_LINKS = $(Opt) -pthread -shared\
			   -Wl,--version-script,"$(MATLAB_ROOT)extern/lib/glnxa64/mexFunction.map" \
			   -Wl,--no-undefined -lblas -llapack
			   

			   
CXX_LIBS = -Wl,--no-undefined -Wl,-rpath-link,"$(MATLAB_ROOT)bin/glnxa64" \
		   -L$(MATLAB_ROOT)bin/glnxa64 -lmx -lmex -lmat -lm -fopenmp 
		   
		   

########################## Mesh ###########################		   
TriangleWrapper    = src/MeshWrapper/TriangleWrapper
TriangleWrapperOut = class/TriangleMesh/private

$(TriangleWrapper)/tTriangleInfo.o: $(TriangleWrapper)/tTriangleInfo.cpp $(TriangleWrapper)/tTriangleInfo.h path
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $< -o $@

$(TriangleWrapperOut)/TriangleWrapper.mexa64: $(TriangleWrapper)/tTriangleInfo.o
	$(CXX) $(MATLAB_LINKS) -o $@ $< $(CXX_LIBS) -ltriangle && rm $(TriangleWrapper)/tTriangleInfo.o
###########################################################	
	
####################### Function Space ####################
FunctionWrapper = src/FunctionSpace
FunctionWrapperOut = class/FunctionSpace/private

$(FunctionWrapper)/FunctionSpace.o: $(FunctionWrapper)/FunctionSpace.cpp $(FunctionWrapper)/FunctionSpace.h path
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $< -o $@
	
$(FunctionWrapperOut)/FunctionSpaceWrapper.mexa64: $(FunctionWrapper)/FunctionSpace.o
	$(CXX) $(MATLAB_LINKS) -o $@ $< $(CXX_LIBS) && rm $(FunctionWrapper)/FunctionSpace.o
###########################################################

########################### Mode ##########################
ModeWrapper = src/ModeWrapper
ModeWrapperOut = class/QuadMode/private

$(ModeWrapper)/ModeWrapper.o: $(ModeWrapper)/tModeInfo.cpp $(ModeWrapper)/tModeInfo.h path
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $< -o $@
	
$(ModeWrapperOut)/ModeWrapper.mexa64: $(ModeWrapper)/ModeWrapper.o
	$(CXX) $(MATLAB_LINKS) -o $@ $< $(CXX_LIBS) -lquadmath && rm  $(ModeWrapper)/ModeWrapper.o

###########################################################

######################## Metis ############################
MetisWrapper = src/MetisWrapper
MetisWrapperOut = utility/MeshPartition

$(MetisWrapper)/MetisWrapper.o: $(MetisWrapper)/metismex.c path
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $< -o $@
	
$(MetisWrapperOut)/MetisPartition.mexa64: $(MetisWrapper)/MetisWrapper.o
	$(CXX) $(MATLAB_LINKS) -o $@ $< $(CXX_LIBS) -lmetis && rm  $(MetisWrapper)/MetisWrapper.o
###########################################################

##################### Form ################################
FormWrapper = src/FormWrapper
FormWrapperOut = class/FormBuilder/private


$(FormWrapper)/FormWrapper.o: $(FormWrapper)/Assembler.cpp path
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $< -o $@
	
$(FormWrapperOut)/FormWrapper.mexa64: $(FormWrapper)/FormWrapper.o
	$(CXX) $(MATLAB_LINKS) -o $@ $< $(CXX_LIBS) && rm  $(FormWrapper)/FormWrapper.o

###########################################################

##################### BC ##################################
BCWrapper = src/BCWrapper
BCWrapperOut = class/BC/private

$(BCWrapper)/BCWrapper.o: $(BCWrapper)/Boundary.cpp path
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $< -o $@
	
$(BCWrapperOut)/BCWrapper.mexa64: $(BCWrapper)/BCWrapper.o
	$(CXX) $(MATLAB_LINKS) -o $@ $< $(CXX_LIBS) && rm  $(BCWrapper)/BCWrapper.o

###########################################################
path:
	mkdir -p ./class/TriangleMesh/private
	mkdir -p ./class/QuadMode/private
	mkdir -p ./class/FunctionSpace/private
	mkdir -p ./class/FormBuilder/private
	mkdir -p ./class/BC/private 
	mkdir -p ./utility/MeshPartition
	
all:$(TriangleWrapperOut)/TriangleWrapper.mexa64 \
$(FunctionWrapperOut)/FunctionSpaceWrapper.mexa64 \
$(ModeWrapperOut)/ModeWrapper.mexa64 \
$(MetisWrapperOut)/MetisPartition.mexa64 \
$(FormWrapperOut)/FormWrapper.mexa64 \
$(BCWrapperOut)/BCWrapper.mexa64



clean:
	rm -f $(TriangleWrapperOut)/TriangleWrapper.mexa64
	rm -f $(FunctionWrapperOut)/FunctionSpaceWrapper.mexa64
	rm -f $(ModeWrapperOut)/ModeWrapper.mexa64
	rm -f $(MetisWrapperOut)/MetisPartition.mexa64
	rm -f $(FormWrapperOut)/FormBuilder.mexa64
	rm -rf $(BCWrapperOut)/BCWrapper.mexa64


