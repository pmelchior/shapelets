#!/usr/bin/env python
import os,sys

# path prefix for installation location
# typical choices are "$(HOME)" or "/usr"
# it is only used for finding the paths when doing a 'make install'
# libshapelens.[a,so] will go to PREFIX/lib
# contents of include will go to PREFIX/include/shapelens
PREFIX = "$(HOME)"

# choose the appropriate values for the compiler and flags
COMPILER = "g++"
FLAGS = "-ansi -g -DNDEBUG -Wno-deprecated -O3 -march=pentium4"

# specify non-standard include locations for headers and libraries
# in particular, give location of numla headers
INCLUDES = ["/usr/include/mysql","$(HOME)/include/numla","$(HOME)/include/atlas","$(HOME)/include"]
LIBS = []

# create shared library
SHARED = True

# set this to False if you don't want an automatically generated documentation
# if you set it to True, doxygen should be installed on your machine.
DOCUMENTATION = True

# set this to False if you do not have FFTW3 installed on the machine,
# otherwise set it to True
FFTW3 = True

# set this to False if you don't want a data base backend compiled,
# set it to "MySQL" if you want to have to MySQL backend compiled
SHAPELETDB = "MySQL"

##########################################
# no changes necessary beyond this point #
##########################################

SUBDIR = os.environ.get('SUBDIR')
MAKEFILE = "Makefile"
if (SUBDIR != None):
	MAKEFILE = MAKEFILE+"."+SUBDIR
else:
	SUBDIR = ""

# insert PATH and CFLAGS in Makefile
print "@@@ Creating new "+ MAKEFILE +" @@@"
infile = open('Makefile.in','r') 
makefile = infile.read()
infile.close()
general = "\nPREFIX = " + PREFIX + "\n"

targets = "all: lib"
doctarget = ""
if (SHARED == True):
        targets = targets + " shared"
if (DOCUMENTATION == True):
	targets = targets + " docs"
	doctarget = "docs: $(HEADERS)\n\tdoxygen Doxyfile"
if (SHAPELETDB == "MySQL"):
	FLAGS = FLAGS + " -DBIG_JOINS=1 -DSHAPELETDB=" + SHAPELETDB
if (FFTW3 == True):
	FLAGS = FLAGS + " -DHAS_FFTW3"

# loop thru paths and append every memeber of INCLUDES to otherincl
otherincl = ""
otherlibs = ""
speciallibs = ""
for i in INCLUDES:
	otherincl = otherincl + " -I"+i
for l in LIBS:
	otherlibs = otherlibs + " -L"+l
if FFTW3:
	speciallibs = speciallibs + " -lfftw3"
if SHAPELETDB == "MySQL":
	speciallibs = speciallibs + " -lmysqlclient"

makefile = makefile.replace("???GENERAL???",general)
makefile = makefile.replace("???COMPILER???",COMPILER)
makefile = makefile.replace("???CFLAGS???",FLAGS)
makefile = makefile.replace("???OTHERINCL???",otherincl)
makefile = makefile.replace("???OTHERLIBS???",otherlibs)
makefile = makefile.replace("???SPECIALLIBS???",speciallibs)
makefile = makefile.replace("???TARGETS???",targets)
makefile = makefile.replace("???DOCTARGET???",doctarget)	

outfile = open(MAKEFILE,'w')
outfile.write(makefile)
outfile.close()

# create directories
if not os.path.isdir("lib/"+SUBDIR):
	print "@@@ Creating directory lib/"+SUBDIR + " @@@"
	os.system("mkdir -p lib/"+SUBDIR)
if DOCUMENTATION == True and not os.path.isdir("docs"):
	print "@@@ Creating directory doc @@@"
	os.system("mkdir -p docs")
