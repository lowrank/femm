CXX = g++
CC  = gcc
FF  = gfortran
Opt = -Ofast

include Makefile.in

CXX_FLAGS = -DMATLAB_MEX_FILE -std=c++11 -fopenmp -march=native \
			-D_GNU_SOURCE -fPIC -fno-omit-frame-pointer -Wno-write-strings -pthread\
			-DMX_COMPAT_32 $(Opt) -DNDEBUG -fopenmp -ffast-math
			
CXX_INCLUDE = -I./include/mexplus \
              -I./include/ \
              -I$(MATLAB_ROOT)extern/include \
              -I$(MATLAB_ROOT)simulink/include
              
MATLAB_LINKS = $(Opt) -pthread -shared\
			   -Wl,--version-script,$(MATLAB_ROOT)extern/lib/glnxa64/mexFunction.map \
			   -Wl,--no-undefined -lblas -llapack
			   
CXX_LIBS = -Wl,-rpath-link,$(MATLAB_ROOT)bin/glnxa64 \
		   -L$(MATLAB_ROOT)bin/glnxa64 -lmx -lmex -lmat -lm -fopenmp
		   
###########################################################		   
TriangleWrapper    = src/MeshWrapper/TriangleWrapper
TriangleWrapperOut = class/TriangleMesh/private

$(TriangleWrapper)/tTriangleInfo.o: $(TriangleWrapper)/tTriangleInfo.cpp $(TriangleWrapper)/tTriangleInfo.h
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $< -o $@

$(TriangleWrapperOut)/TriangleWrapper.mexa64: $(TriangleWrapper)/tTriangleInfo.o
	$(CXX) $(MATLAB_LINKS) -o $@ $< $(CXX_LIBS) -ltriangle && rm $(TriangleWrapper)/tTriangleInfo.o
###########################################################	
	
###########################################################
FunctionWrapper = src/FunctionSpace
FunctionWrapperOut = class/FunctionSpace/private

$(FunctionWrapper)/FunctionSpace.o: $(FunctionWrapper)/FunctionSpace.cpp $(FunctionWrapper)/FunctionSpace.h
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $< -o $@
	
$(FunctionWrapperOut)/FunctionSpaceWrapper.mexa64: $(FunctionWrapper)/FunctionSpace.o
	$(CXX) $(MATLAB_LINKS) -o $@ $< $(CXX_LIBS) && rm $(FunctionWrapper)/FunctionSpace.o
###########################################################

	
all:$(TriangleWrapperOut)/TriangleWrapper.mexa64 $(FunctionWrapperOut)/FunctionSpaceWrapper.mexa64