####################################################
# Copyright 2010,2011,2012 Alois Schloegl 
# This is part of the NaN-toolbox - a statistics and machine learning toolbox for data with and without missing values. 
# http://pub.ist.ac.at/~schloegl/matlab/NaN/
####################################################

### modify directories according to your needs 

# Octave - global install  (e.g. from debian package)
# OCTAVE_VERSION = 

# Matlab configuration
MATLABDIR  = /usr/local/MATLAB/R2010b
# comment the following line if you use MATLAB on 32-bit operating system
MEX_OPTION += -largeArrayDims

# Mingw crosscompiler: available at http://www.nongnu.org/mingw-cross-env/
CROSS   = $(HOME)/src/mxe/usr/bin/i686-pc-mingw32
CROSS64 = $(HOME)/src/mce-w64/usr/bin/x86_64-static-mingw32
# include directory for Win32-Matlab include
W32MAT_INC = -I$(HOME)/bin/win32/Matlab/R2010b/extern/include/
W64MAT_INC = -I$(HOME)/bin/win64/Matlab/R2010b/extern/include/
# path to GNUMEX libraries, available from here http://sourceforge.net/projects/gnumex/
GNUMEX   = $(HOME)/bin/win32/gnumex
GNUMEX64 = $(HOME)/bin/win64/gnumex
# building gnumex64 was difficult, these hints were quite useful: 
# http://sourceforge.net/mailarchive/forum.php?thread_name=AANLkTinZvxgC9ezp2P3UCX_a7TAUYuVsp2U40MQUV6qr%40mail.gmail.com&forum_name=gnumex-users
# Instead of building "mex shortpath.c" and "mex uigetpath.c", I used empty m-functions within argout=argin; 
####################################################
W32MAT_INC += -I$(subst /usr/bin/,/usr/,$(CROSS))/include/
W64MAT_INC += -I$(subst /usr/bin/,/usr/,$(CROSS64))/include/
LDLIBS_W32 = $(subst /usr/bin/,/usr/,$(CROSS))/lib/
LDLIBS_W64 = $(subst /usr/bin/,/usr/,$(CROSS64))/lib/
W32_LIBS  = $(LDLIBS_W32)liblapack.a
W64_LIBS  = $(LDLIBS_W64)liblapack.a
W32_LIBS += $(LDLIBS_W32)libblas.a
W64_LIBS += $(LDLIBS_W64)libblas.a



CC  	= gcc
CXX 	= g++
CFLAGS 	= -fopenmp -Wall -Wextra -Wconversion -O2 -fPIC 
OCTMEX	= mkoctfile$(OCTAVE_VERSION) --mex
RM      = rm

MEX_OPTION  += -lgomp CC\#$(CXX) CXX\#$(CXX) CFLAGS\#"$(CFLAGS)" CXXFLAGS\#"$(CFLAGS)"
MATMEX      = $(MATLABDIR)/bin/mex $(MEX_OPTION)

PROGS = covm_mex.mex sumskipnan_mex.mex 

mex4o octave: $(PROGS)

ifneq ( , $(ls -1 $(MATLABDIR) ) )
MEX_EXT = $(shell $(MATLABDIR)/bin/mexext)
mex4m matlab: $(patsubst %.mex, %.$(MEX_EXT), $(PROGS))
endif
mexw32 win32:  $(patsubst %.mex, %.mexw32, $(PROGS))
mexw64 win64:  $(patsubst %.mex, %.mexw64, $(PROGS))
all:    octave win32 win64 mex4m

clean:
	-$(RM) *.o *.obj *.o64 core octave-core *.oct *~ *.mex* 


#########################################################

#	Octave, MATLAB on Linux 
#########################################################
%.oct: %.cc
	mkoctfile$(OCTAVE_VERSION) "$<"

%.mex: %.cpp
	$(OCTMEX) "$<" -lgomp -llapack -lblas     ## Octave
%.$(MEX_EXT): %.cpp
	$(MATMEX) "$<" -lgomp -llapack -lblas     ## Matlab


#########################################################
#	MATLAB/WIN32
#########################################################
%.obj: %.cpp
	$(CROSS)-$(CXX) -fopenmp -c -DMATLAB_MEX_FILE -x c++ -o "$@" $(W32MAT_INC) -O2 -DMX_COMPAT_32 "$<"
%.obj: %.c
	$(CROSS)-$(CXX) -fopenmp -c -DMATLAB_MEX_FILE -x c++ -o "$@" $(W32MAT_INC) -O2 -DMX_COMPAT_32 "$<" 

%.mexw32: %.obj
	$(CROSS)-$(CXX) -shared $(GNUMEX)/mex.def -o "$@" -L$(GNUMEX) -s "$<" -llibmx -llibmex -llibmat -lcholmod -lgomp -lpthread -L$(LDLIBS_W32) -lblas -llapack


#########################################################
#	MATLAB/WIN64
#########################################################

## ToDO: fix OpenMP support: currently -fopenmp causes Matlab to crash
%.o64: %.cpp
	$(CROSS64)-$(CXX) -c -DMATLAB_MEX_FILE -x c++ -o "$@" $(W64MAT_INC) -O2 "$<"
%.o64: %.c
	$(CROSS64)-$(CXX) -c -DMATLAB_MEX_FILE -x c++ -o "$@" $(W64MAT_INC) -O2 "$<" 

%.mexw64: %.o64
	$(CROSS64)-$(CXX) -shared $(GNUMEX64)/mex.def -o "$@" -L$(GNUMEX64) -s "$<" -llibmx -llibmex -llibmat -lcholmod -lgomp -lpthread -L$(LDLIBS_W64) -lblas -llapack
	
	

