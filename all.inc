#General make file settings
CC = cc
CPP = c++
LIBTOOL = ar -r

#platform-specific options
#Compaq compiler C++ options:
#PCPPOPTS=-D__USE_STD_IOSTREAM -tweak

#MPICH-2 specific option
#PCPPOPTS=-DMPICH_IGNORE_CXX_SEEK

#Altivec specific C options:
#PCOPTS=-maltivec -D _OPT_ALTIVEC_

#Pthread linking support options
PTHREADLD=-lpthread

#OpenMP specific compiler options
#OPENMPOPT=-fopenmp

#OpenMP specific linker options
#OPENMPLD=-fopenmp

#optimization options
OPTIMOPT=-g
#OPTIMOPT=-O

#Poll option
#POLLOPT=-D__USE_POLL__

#Static linking option (if needed)
#STATIC=-static

#C++ 11 support on
STD_OPT = -std=c++14


# Auxilary headers
INITTESTS=$(ROOT)gtest/include
SNOWGOOSE=$(ROOT)../snowgoose
COMPI=$(ROOT)../COMPI
OPTITEST=$(ROOT)../OPTITEST
LOCSEARCH=$(ROOT)../LOCSEARCH
LATTICE=$(ROOT)../2DLATTICE
CRYSTAL=$(ROOT)../CRYSTAL
NUC=$(ROOT)../NUC
EXTRA_INC = -I$(INITTESTS) -I$(SNOWGOOSE) -I$(COMPI) -I$(OPTITEST) -I$(LOCSEARCH) -I$(LATTICE) -I$(CRYSTAL) -I$(NUC)


#common options
COMOPTS = $(STD_OPT)\
	   -I$(ROOT)\
           $(EXTRA_INC)\
       $(OPENMPOPT)\
       $(OPTIMOPT)\
       $(POLLOPT)
#C options
COPTS=$(COMOPTS) $(ADDIT_INCS) $(PCOPTS)

#C++ options
CPPOPTS= $(COMOPTS) $(ADDIT_INCS) $(PCPPOPTS)

#Linker options for static linking
LDFLAGS = -lm -ldl $(OPENMPLD) $(PTHREADLD) $(STATIC) -L$(ROOT)gtest/lib -lgtest


MFLAGS = MPICC=$(MPICC) COPTS="$(COPTS)" CPPOPTS="$(CPPOPTS) LIBTOOL="$(LIBTOOL)"

GTEST_LIB = gtest/lib
GTEST_JSON_PATH= "$(ROOT)../OPTITEST/testfuncs/funcdesc.json"

all: $(BINS) $(LIBS)

tests:
	@for i in $(TESTS); do ./$$i $(GTEST_JSON_PATH); continue; done

# Uncomment this if you want testing finishes after first finished test
#	@for i in $(TESTS); do if ./$$i > /dev/null; then echo TEST PASSED; continue; else echo TEST FAILED; return 1; fi done
	

$(BINS): $(LIBS) $(EXTLIBS)

$(LIBS): $(OBJS)
	$(LIBTOOL) $@ $(OBJS)

indent: 
	indent *.cpp ; indent *.hpp
clean: 
	rm -f *.exe *.o *.a *~ *.log deps.inc

dep:
	$(CC) $(COPTS) -MM -c *.c > deps.inc;\
	$(CPP) $(CPPOPTS) -MM -c *.cpp >> deps.inc;\
        true

.o.exe:
	$(CPP) $(OPTS) -o $@ $< $(LIBS) $(EXTLIBS) $(LDFLAGS)

.cpp.o:
	$(CPP) $(CPPOPTS) -c $<

.c.o:
	$(CC) $(COPTS) -c $<

.SUFFIXES:
.SUFFIXES: .o .a .cpp .c .exe
