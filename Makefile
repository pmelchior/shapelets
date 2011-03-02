INCLPATH = include
SRCPATH = src
LIBPATH = lib
DOCPATH = doc
PROGSRCPATH = progs
LIBNAME = shapelens

SHAPELETSSRCPATH = $(SRCPATH)/shapelets
SHAPELETSINCLPATH = $(INCLPATH)/shapelets
SHAPELETSSRC = $(wildcard $(SHAPELETSSRCPATH)/*.cc)
SHAPELETSOBJECTS = $(SHAPELETSSRC:$(SHAPELETSSRCPATH)/%.cc=$(LIBPATH)/%.o)

LENSINGSRCPATH = $(SRCPATH)/lensing
LENSINGINCLPATH = $(INCLPATH)/lensing
LENSINGSRC = $(wildcard $(LENSINGSRCPATH)/*.cc)
LENSINGOBJECTS = $(LENSINGSRC:$(LENSINGSRCPATH)/%.cc=$(LIBPATH)/%.o)

FRAMESRCPATH = $(SRCPATH)/frame
FRAMEINCLPATH = $(INCLPATH)/frame
FRAMESRC = $(wildcard $(FRAMESRCPATH)/*.cc)
FRAMEOBJECTS = $(FRAMESRC:$(FRAMESRCPATH)/%.cc=$(LIBPATH)/%.o)

COMMONSRC = $(wildcard $(SRCPATH)/*.cc)
COMMONOBJECTS = $(COMMONSRC:$(SRCPATH)/%.cc=$(LIBPATH)/%.o)

UTILSSRCPATH = $(SRCPATH)/utils
UTILSINCLPATH = $(INCLPATH)/utils
UTILSSRC = $(wildcard $(SRCPATH)/utils/*.cc)
UTILSOBJECTS = $(UTILSSRC:$(SRCPATH)/utils/%.cc=$(LIBPATH)/%.o)

ITALIBSLIBPATH = $(ITALIBSPATH)/lib
PROGPATH = $(ITALIBSPATH)/bin
PROGS = $(wildcard $(PROGSRCPATH)/*.cc)
PROGSOBJECTS = $(PROGS:$(PROGSRCPATH)/%.cc=$(PROGPATH)/%)

HEADERS = $(wildcard $(INCLPATH)/*.h) $(wildcard $(SHAPELETSINCLPATH)/*.h) $(wildcard $(FRAMEINCLPATH)/*.h) $(wildcard $(LENSINGINCLPATH)/*.h) $(wildcard $(UTILSINCLPATH)/*.h)

# which OS
UNAME := $(shell uname)

CC = g++
#compilation flags
CFLAGS = -ansi -g $(SPECIALFLAGS) -DBIG_JOINS=1
ifneq ($(UNAME),Linux)
	CFLAGS += -bind_at_load
endif

ifneq (,$(findstring HAS_ATLAS,$(SPECIALFLAGS)))
	LIBS = -lshapelens -lgsl -llapack_atlas -latlas -llapack -lcfitsio -lfftw3 -lsqlite3
else
	LIBS = -lshapelens -lgsl -lgslcblas -lcfitsio -lfftw3 -lsqlite3
endif

AR = ar
ARFLAGS = -sr
ifeq ($(UNAME),Linux)
	SHAREDFLAGS = -shared -fPIC 
	LIBEXT = so
else
	SHAREDFLAGS = -dynamiclib -fPIC
	LIBEXT = dylib
endif

all: $(LIBPATH) $(DOCPATH) library shared

$(LIBPATH):
	mkdir -p $(LIBPATH)

$(DOCPATH):
	mkdir -p $(DOCPATH)

.PHONY: clean

clean: 
	rm -f $(SHAPELETSOBJECTS) $(FRAMEOBJECTS) $(LENSINGOBJECTS) $(MODELFITOBJECTS) $(COMMONOBJECTS) $(UTILSOBJECTS)
	rm -f $(LIBPATH)/lib$(LIBNAME).a
	rm -f $(LIBPATH)/lib$(LIBNAME).$(LIBEXT)	
	rm -f $(DOCPATH)/*

cleanshared:
	rm -f $(LIBPATH)/lib$(LIBNAME).$(LIBEXT)

cleanshapelets:
	rm -f $(SHAPELETSOBJECTS)
	rm -f $(LIBPATH)/lib$(LIBNAME).a
	rm -f $(LIBPATH)/lib$(LIBNAME).$(LIBEXT)

cleanlensing:
	rm -f $(LENSINGOBJECTS)
	rm -f $(LIBPATH)/lib$(LIBNAME).a
	rm -f $(LIBPATH)/lib$(LIBNAME).$(LIBEXT)

cleanframe:
	rm -f $(FRAMEOBJECTS)
	rm -f $(LIBPATH)/lib$(LIBNAME).a
	rm -f $(LIBPATH)/lib$(LIBNAME).$(LIBEXT)

cleancommon:
	rm -f $(COMMONOBJECTS)
	rm -f $(LIBPATH)/lib$(LIBNAME).a
	rm -f $(LIBPATH)/lib$(LIBNAME).$(LIBEXT)

cleanprogs:
	rm -f $(PROGSOBJECTS)

library: $(LIBPATH)/lib$(LIBNAME).a

shared: $(LIBPATH)/lib$(LIBNAME).$(LIBEXT)

install: library shared
	mkdir -p $(ITALIBSLIBPATH)
	cp $(LIBPATH)/lib$(LIBNAME).a $(ITALIBSLIBPATH)
	cp $(LIBPATH)/lib$(LIBNAME).$(LIBEXT) $(ITALIBSLIBPATH)
	mkdir  -p $(ITALIBSPATH)/include/$(LIBNAME)
	cd $(INCLPATH) && find . -type f -name '*.h' -exec  cp --parents {} $(ITALIBSPATH)/include/$(LIBNAME)/ \; && cd ../
	mkdir -p $(PROGPATH)

progs: $(PROGSOBJECTS)

docs: $(HEADERS)
	doxygen Doxyfile

$(LIBPATH)/lib$(LIBNAME).a: $(SHAPELETSOBJECTS) $(FRAMEOBJECTS) $(LENSINGOBJECTS) $(COMMONOBJECTS) $(UTILSOBJECTS)
	$(AR) $(ARFLAGS) $@ $?

ifeq ($(UNAME),Linux)
$(LIBPATH)/lib$(LIBNAME).$(LIBEXT): $(SHAPELETSOBJECTS) $(FRAMEOBJECTS) $(LENSINGOBJECTS) $(COMMONOBJECTS) $(UTILSOBJECTS)
	rm -f $(LIBPATH)/lib$(LIBNAME).$(LIBEXT)	
	$(CC) $(SHAREDFLAGS) -o $@ $^
else
$(LIBPATH)/lib$(LIBNAME).$(LIBEXT): $(SHAPELETSOBJECTS) $(FRAMEOBJECTS) $(LENSINGOBJECTS) $(COMMONOBJECTS) $(UTILSOBJECTS)
	rm -f $(LIBPATH)/lib$(LIBNAME).$(LIBEXT)	
	$(CC) $(SHAREDFLAGS) $(CFLAG_LIBS) -o $@ $^ $(LIBS)
endif

$(LIBPATH)/%.o: $(SHAPELETSSRCPATH)/%.cc
	$(CC) $(CFLAGS) -c $< -o $@

$(LIBPATH)/%.o: $(FRAMESRCPATH)/%.cc
	$(CC) $(CFLAGS) -c $< -o $@

$(LIBPATH)/%.o: $(LENSINGSRCPATH)/%.cc
	$(CC) $(CFLAGS) -c $< -o $@

$(LIBPATH)/%.o: $(UTILSSRCPATH)/%.cc
	$(CC) $(CFLAGS) -c $< -o $@

$(LIBPATH)/%.o: $(SRCPATH)/%.cc
	$(CC) $(CFLAGS) -c $< -o $@

$(PROGPATH)/%: $(PROGSRCPATH)/%.cc
	$(CC) $(CFLAGS) $(CFLAG_LIBS) $< -o $@ $(LIBS)
