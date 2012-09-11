# --- variables defined by configure ---

SRC = ./src
PREFIX = x86_64-Darwin
LIBDIRSUFFIX = 

EXE = 
DEF = -D
NOUNDERSCORE = 1

FC = gfortran
FFLAGS = -ff2c -O0 -g -ffixed-line-length-none     -fno-range-check -ff2c  $(DEF)QUAD=0 $(DEF)U77EXT=0

CC = clang
CFLAGS = -O3 -g -fomit-frame-pointer -ffast-math -Wall -m32   -DNOUNDERSCORE=$(NOUNDERSCORE) -DBIGENDIAN=0

CXX = clang
CXXFLAGS = -O3 -g -fomit-frame-pointer -ffast-math -Wall -m32  

ML = 0
MCC = mcc
MCFLAGS =   -b32

LDFLAGS = -lpthread -L/usr/local/Cellar/gfortran/4.2.4-5666.3/bin/../lib/gcc/i686-apple-darwin11/4.2.1/x86_64 -L/usr/lib/gcc/i686-apple-darwin11/4.2.1/x86_64 -L/usr/lib/gcc/i686-apple-darwin11/4.2.1/x86_64 -L/usr/lib/i686-apple-darwin11/4.2.1 -L/usr/local/Cellar/gfortran/4.2.4-5666.3/bin/../lib/gcc/i686-apple-darwin11/4.2.1 -L/usr/local/Cellar/gfortran/4.2.4-5666.3/bin/../lib/gcc -L/usr/lib/gcc/i686-apple-darwin11/4.2.1 -L/usr/lib/gcc/i686-apple-darwin11/4.2.1 -L/usr/local/Cellar/gfortran/4.2.4-5666.3/bin/../lib/gcc/i686-apple-darwin11/4.2.1/../../.. -L/usr/lib/gcc/i686-apple-darwin11/4.2.1/../../../i686-apple-darwin11/4.2.1 -L/usr/lib/gcc/i686-apple-darwin11/4.2.1/../../.. -lgfortranbegin -lgfortran -lSystem -lSystem -m32

AR = ar
RANLIB = ranlib
NM = nm
DLLTOOL = 

PARALLEL = -j 2

# --- end defs by configure ---


BLD = build$(QUADSUFFIX)

LIBDIR = $(PREFIX)/lib$(LIBDIRSUFFIX)
INCLUDEDIR = $(PREFIX)/include
BINDIR = $(PREFIX)/bin

LIB = libooptools$(QUADSUFFIX).a
FE = lt$(QUADSUFFIX)$(EXE)
MFE = LoopTools$(QUADSUFFIX)$(EXE)
INCLUDE = $(BLD)/looptools.h $(BLD)/clooptools.h

ARGS = $(PARALLEL) \
  LIB="$(LIB)" \
  FE="$(FE)" \
  MFE="$(MFE)" \
  EXE="$(EXE)" \
  DEF="$(DEF)" \
  NOUNDERSCORE="$(NOUNDERSCORE)" \
  XFC="$(FC$(QUADSUFFIX)) $(FFLAGS) $(FFLAGS-quad) -I." \
  F90="$(F90)" \
  CC="$(CC)" \
  CFLAGS="$(CFLAGS) $(CFLAGS-quad)" \
  CXX="$(CXX)" \
  CXXFLAGS="$(CXXFLAGS)" \
  ML="$(ML)" \
  MCC="$(MCC)" \
  MCFLAGS="$(MCFLAGS)" \
  AR="$(AR)" \
  RANLIB="$(RANLIB)" \
  NM="$(NM)" \
  DLLTOOL="$(DLLTOOL)" \
  LDFLAGS="$(LDFLAGS)" \
  LIBPATH="$(LIBPATH)" \
  OBJS-quad="$(OBJS-quad)"


default all lib frontend mma: force
	cd $(BLD) && $(MAKE) $(ARGS) $@

install: lib frontend
	-mkdir $(PREFIX)
	-mkdir $(LIBDIR) $(BINDIR) $(INCLUDEDIR)
	cp -p $(BLD)/$(LIB) $(LIBDIR)
	cp -p $(INCLUDE$(QUADSUFFIX)) $(INCLUDEDIR)
	strip $(BLD)/$(FE)
	cp -p $(BLD)/fcc $(BLD)/$(FE) $(BINDIR)
	test ! -f $(BLD)/$(MFE) || { strip $(BLD)/$(MFE) ; cp -p $(BLD)/$(MFE) $(BINDIR); }

force: $(BLD)/timestamp

$(BLD)/timestamp:
	-mkdir $(BLD)
	find "`cd $(SRC) ; pwd`" -follow -exec ln -sf {} $(BLD) \;
	touch $(BLD)/timestamp

clean:
	rm -fr $(BLD)

