#!/usr/bin/python
import os

# change this according to your machine and paths
# compiled lib and progs will go into lib/$SUBDIR and progs/$SUBDIR
INSTALLPATH = "."
SUBDIR = "i686/pentium4"	
SHAPELENSPATH = "$(HOME)/code/shapelens"
NUMLAINCLPATH = "$(HOME)/code/numla"
OTHERINCL = ["$(HOME)/code/boost-local","$(HOME)/include","$(HOME)/code/ATLAS/include"]
OTHERLIB = ["$(HOME)/lib"]

# choose the appropriate values for the compiler and flags
COMPILER = "g++"
FLAGS = "-ansi -g -DNDEBUG -Wno-deprecated -O3 -march=pentium4"

# create shared library
SHARED = True

# set this to False if you don't want an automatically generated documentation
# if set to True, doxygen should be installed on your machine.
DOCUMENTATION = True



##########################################
# no changes necessary beyond this point #
##########################################

# insert PATH and CFLAGS in Makefile
infile = open('Makefile.in','r') 
makefile = infile.read()
infile.close()
paths = "INSTALLPATH = "+INSTALLPATH+"\nSHAPELENSPATH = "+SHAPELENSPATH+"\nNUMLAINCLPATH = "+NUMLAINCLPATH

targets = "all: lib"
doctarget = ""
if (SHARED == True):
        targets = targets + " shared"
if (DOCUMENTATION == True):
	targets = targets + " docs"
	doctarget = "docs: $(HEADERS)\n\tdoxygen Doxyfile"
targets = targets + " progs"

otherincl = ""
for i in range(len(OTHERINCL)):
	otherincl = otherincl + " -I"+OTHERINCL[i]
otherlib = ""
for i in range(len(OTHERLIB)):
        otherlib = otherlib + " -L"+OTHERLIB[i]

makefile = makefile.replace("???PATHS???",paths)
makefile = makefile.replace("???SUBDIR???",SUBDIR)
makefile = makefile.replace("???COMPILER???",COMPILER)
makefile = makefile.replace("???CFLAGS???",FLAGS)
makefile = makefile.replace("???OTHERINCL???",otherincl)
makefile = makefile.replace("???OTHERLIB???",otherlib)
makefile = makefile.replace("???TARGETS???",targets)
makefile = makefile.replace("???DOCTARGET???",doctarget)	

outfile = open('Makefile','w')
outfile.write(makefile)
outfile.close()


# insert INPUT in Doxyfile
if (DOCUMENTATION == True):
	infile = open('Doxyfile.in','r')
	doxyfile = infile.read()
	infile.close()

	input = NUMLAINCLPATH+" "+SHAPELENSPATH+"/include"
	doxyfile = doxyfile.replace("???INPUT???",input)

	outfile = open('Doxyfile','w')
	outfile.write(doxyfile)
	outfile.close()

# create directories
os.system("mkdir -p progs/"+SUBDIR)
os.system("mkdir -p lib/"+SUBDIR)
if (DOCUMENTATION == True):
	os.system("mkdir -p docs")
