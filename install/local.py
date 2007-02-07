#!/usr/bin/python
import os

# change this according to you machine and paths
INSTALLPATH = "."
SRCPATH = "../shapelens++"
NUMLAINCLDIR = "../numla"
ATLASINCLDIR = "../ATLAS/include"
MARCH="pentium4"
OPTIMIZE="3"

# no changes necessary beyond this point
# insert PATH and CFLAGS in Makefile
infile = open('Makefile.in','r') 
makefile = infile.read()
infile.close()
paths = "INSTALLPATH = "+INSTALLPATH+"\nSRCPATH = "+SRCPATH+"\nNUMLAINCLDIR = "+NUMLAINCLDIR+"\nATLASINCLDIR = "+ATLASINCLDIR

makefile = makefile.replace("???PATHS???",paths)
makefile = makefile.replace("???OPTIMIZE???",OPTIMIZE)
makefile = makefile.replace("???MARCH???",MARCH)

outfile = open('Makefile','w')
outfile.write(makefile)
outfile.close()


# insert INPUT in Doxyfile
infile = open('Doxyfile.in','r')
doxyfile = infile.read()
infile.close()

input = NUMLAINCLDIR+" "+SRCPATH+"/shapelets/include "+SRCPATH+"/lensing/include "+SRCPATH+"/frame/include"
doxyfile = doxyfile.replace("???INPUT???",input)

outfile = open('Doxyfile','w')
outfile.write(doxyfile)
outfile.close()

# create lib and progs directories
arch = os.popen("uname -m").readline()
os.system("mkdir -p lib/"+str(arch))
os.system("mkdir -p progs/"+str(arch))
