#!/usr/bin/python
import os

# change this according to you machine and paths
INSTALLPATH = "."
SRCPATH = "../shapelens"
NUMLAINCLDIR = "../numla"
OTHERINCL = ["$(HOME)/include","../ATLAS/include"]
OTHERLIB = ["$(HOME)/lib","../ATLAS/lib/Linux_P4SSE3_2"]

MARCH="pentium4"
OPTIMIZE="3"

# set this to False if you don't want a automatically generated documentation
# if set to True, doxygen should be installed on your machine.
DOCUMENTATION = True

# no changes necessary beyond this point
# insert PATH and CFLAGS in Makefile
infile = open('Makefile.in','r') 
makefile = infile.read()
infile.close()
paths = "INSTALLPATH = "+INSTALLPATH+"\nSRCPATH = "+SRCPATH+"\nNUMLAINCLDIR = "+NUMLAINCLDIR

targets = "all: shapelets lensing frame progs docs"
doctarget = "docs: $(HEADERS)\n\tdoxygen Doxyfile"
if (DOCUMENTATION==False):
	targets = "all: shapelets lensing frame progs"
	doctarget = ""

otherincl = ""
for i in range(len(OTHERINCL)):
	otherincl = otherincl + " -I"+OTHERINCL[i]
otherlib = ""
for i in range(len(OTHERLIB)):
        otherlib = otherlib + " -L"+OTHERLIB[i]

makefile = makefile.replace("???PATHS???",paths)
makefile = makefile.replace("???OPTIMIZE???",OPTIMIZE)
makefile = makefile.replace("???MARCH???",MARCH)
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

	input = NUMLAINCLDIR+" "+SRCPATH+"/shapelets/include "+SRCPATH+"/lensing/include "+SRCPATH+"/frame/include"
	doxyfile = doxyfile.replace("???INPUT???",input)

	outfile = open('Doxyfile','w')
	outfile.write(doxyfile)
	outfile.close()

# create directories
arch = os.popen("uname -m").read().strip()
os.system("mkdir -p progs/"+arch)
os.system("mkdir -p lib/"+arch)
if (DOCUMENTATION == True):
	os.system("mkdir docs")
