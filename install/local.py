#!/usr/bin/python

# change this according to you machine and paths
INSTALLPATH = "."
SRCPATH = "../shapelens"
NUMLAINCLDIR = "../numla"
OTHERINCL = ["$(HOME)/include","../ATLAS/include"]
OTHERLIB = ["$(HOME)/lib","../ATLAS/lib/Linux_P4SSE3_2"]

MARCH="pentium4"
OPTIMIZE="3"

# no changes necessary beyond this point
# insert PATH and CFLAGS in Makefile
infile = open('Makefile.in','r') 
makefile = infile.read()
infile.close()
paths = "INSTALLPATH = "+INSTALLPATH+"\nSRCPATH = "+SRCPATH+"\nNUMLAINCLDIR = "+NUMLAINCLDIR

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

