INCLPATH = ./include
SRCPATH = ./src
LIBPATH = ./lib/$(SUBDIR)
PROGSRCPATH = ./progs
LIBNAME = shapelens

NUMLAPATH = $(ITALIBSPATH)/include/numla
ATLASLIBPATH = $(ITALIBSPATH)/lib/$(SUBDIR)

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

MODELFITSRCPATH = $(SRCPATH)/modelfit
MODELFITINCLPATH = $(INCLPATH)/modelfit
MODELFITSRC = $(wildcard $(MODELFITSRCPATH)/*.cc)
MODELFITOBJECTS = $(MODELFITSRC:$(MODELFITSRCPATH)/%.cc=$(LIBPATH)/%.o)

COMMONSRC = $(wildcard $(SRCPATH)/*.cc)
COMMONOBJECTS = $(COMMONSRC:$(SRCPATH)/%.cc=$(LIBPATH)/%.o)

UTILSSRCPATH = $(SRCPATH)/utils
UTILSINCLPATH = $(INCLPATH)/utils
UTILSSRC = $(wildcard $(SRCPATH)/utils/*.cc)
UTILSOBJECTS = $(UTILSSRC:$(SRCPATH)/utils/%.cc=$(LIBPATH)/%.o)

ITALIBSLIBPATH = $(ITALIBSPATH)/lib/$(SUBDIR)
PROGPATH = $(ITALIBSPATH)/bin/$(SUBDIR)
PROGS = $(wildcard $(PROGSRCPATH)/*.cc)
PROGSOBJECTS = $(PROGS:$(PROGSRCPATH)/%.cc=$(PROGPATH)/%)

HEADERS = $(wildcard $(INCLPATH)/*.h) $(wildcard $(SHAPELETSINCLPATH)/*.h) $(wildcard $(FRAMEINCLPATH)/*.h) $(wildcard $(MODELFITINCLPATH)/*.h) $(wildcard $(UTILSINCLPATH)/*.h)

CC = g++
CFLAGS = -ansi -g $(SPECIALFLAGS) -DBIG_JOINS=1 -DSHAPELETDB=MySQL -DHAS_FFTW3 -I$(INCLPATH)  -I/usr/include/mysql -I$(NUMLAPATH)
CFLAG_LIBS = -I$(HOME)/include -L$(ATLASLIBPATH) -L$(LIBPATH)
LIBS = -lshapelens -lgsl -lcblas -llapack_atlas -latlas -llapack -lg2c -lcfitsio  -lfftw3 -lmysqlclient

AR = ar
ARFLAGS = -sr
SHAREDFLAGS = -shared -fPIC

.DEFAULT: all

.PHONY: clean

all: lib shared

clean:
	rm -f $(LIBPATH)/*

cleanshared:
	rm -f $(LIBPATH)/lib$(LIBNAME).so

cleanshapelets:
	rm -f $(SHAPELETSOBJECTS)
	rm -f $(LIBPATH)/lib$(LIBNAME).a
	rm -f $(LIBPATH)/lib$(LIBNAME).so

cleanlensing:
	rm -f $(LENSINGOBJECTS)
	rm -f $(LIBPATH)/lib$(LIBNAME).a
	rm -f $(LIBPATH)/lib$(LIBNAME).so

cleanframe:
	rm -f $(FRAMEOBJECTS)
	rm -f $(LIBPATH)/lib$(LIBNAME).a
	rm -f $(LIBPATH)/lib$(LIBNAME).so

cleanmodelfit:
	rm -f $(MODELFITOBJECTS)
	rm -f $(LIBPATH)/lib$(LIBNAME).a
	rm -f $(LIBPATH)/lib$(LIBNAME).so

cleancommon:
	rm -f $(COMMONOBJECTS)
	rm -f $(LIBPATH)/lib$(LIBNAME).a
	rm -f $(LIBPATH)/lib$(LIBNAME).so

cleanprogs:
	rm -f $(PROGSOBJECTS)

lib: $(LIBPATH)/lib$(LIBNAME).a

shared: $(LIBPATH)/lib$(LIBNAME).so

install: lib shared
	mkdir -p $(ITALIBSLIBPATH)
	cp $(LIBPATH)/lib$(LIBNAME).a $(ITALIBSLIBPATH)
	cp $(LIBPATH)/lib$(LIBNAME).so $(ITALIBSLIBPATH)
	mkdir  -p $(ITALIBSPATH)/include/$(LIBNAME)
	cd $(INCLPATH) && find . -type f -name '*.h' -exec  cp --parents {} $(ITALIBSPATH)/include/$(LIBNAME)/ \; && cd ../
	mkdir -p $(PROGPATH)

progs: $(PROGSOBJECTS)

docs: $(HEADERS)
	doxygen Doxyfile

$(LIBPATH)/lib$(LIBNAME).a: $(SHAPELETSOBJECTS) $(FRAMEOBJECTS) $(LENSINGOBJECTS) $(MODELFITOBJECTS) $(COMMONOBJECTS) $(UTILSOBJECTS)
	$(AR) $(ARFLAGS) $@ $?

$(LIBPATH)/lib$(LIBNAME).so: $(SHAPELETSOBJECTS) $(FRAMEOBJECTS) $(LENSINGOBJECTS) $(MODELFITOBJECTS) $(COMMONOBJECTS) $(UTILSOBJECTS)
	$(CC) $(SHAREDFLAGS) -o $@ $?

$(LIBPATH)/%.o: $(SHAPELETSSRCPATH)/%.cc
	$(CC) $(CFLAGS) -c $< -o $@

$(LIBPATH)/%.o: $(FRAMESRCPATH)/%.cc
	$(CC) $(CFLAGS) -c $< -o $@

$(LIBPATH)/%.o: $(MODELFITSRCPATH)/%.cc
	$(CC) $(CFLAGS) -c $< -o $@

$(LIBPATH)/%.o: $(LENSINGSRCPATH)/%.cc
	$(CC) $(CFLAGS) -c $< -o $@

$(LIBPATH)/%.o: $(UTILSSRCPATH)/%.cc
	$(CC) $(CFLAGS) -c $< -o $@

$(LIBPATH)/%.o: $(SRCPATH)/%.cc
	$(CC) $(CFLAGS) -c $< -o $@

$(PROGPATH)/%: $(PROGSRCPATH)/%.cc
	$(CC) $(CFLAGS) $(CFLAG_LIBS) $< -o $@ $(LIBS)
