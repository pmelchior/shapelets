#!/usr/bin/env python
import os,sys

# change this according to your machine and paths
# compiled lib and progs will go into lib/$SUBDIR and progs/$SUBDIR
INSTALLPATH = "."
SUBDIR = "i686/develop"	

# give the name of the file which lists the header and library directories
# if this file does not exist, it will be created
CONFIGFILE = "PATHS"

# choose the appropriate values for the compiler and flags
COMPILER = "g++"
FLAGS = "-ansi -g -DNDEBUG -Wno-deprecated -O3 -march=pentium4"
LIBS = "-lgsl -llapack_atlas -latlas -llapack -lg2c -lcfitsio"

# create shared library
SHARED = True

# set this to False if you don't want an automatically generated documentation
# if set to True, doxygen should be installed on your machine.
DOCUMENTATION = True


##########################################
# no changes necessary beyond this point #
##########################################

# return copy of inlist where every entry is unique
def unique(inlist, keepstr=True):
	typ = type(inlist)
	if not typ == list:
		inlist = list(inlist)
	i = 0
	while i < len(inlist):
		try:
			del inlist[inlist.index(inlist[i], i + 1)]
		except:
			i += 1
	if not typ in (str, unicode):
		inlist = typ(inlist)
	else:
		if keepstr:
			inlist = ''.join(inlist)
	return inlist

# return directory part of path, cutting cutdirs level from the right
def pathTruncated(filename,cutdirs=0):
	path = "/".join(filename.strip().split("/")[0:-1-cutdirs])
	# replace $HOME part of path by "$(HOME)
	home = os.environ["HOME"]
	path = path.replace(home,"$(HOME)")
	return path

# call commandline utility locate to find the path pattern(s)
# if there are multiple directories mathing, it returns the first one
def checkWithLocate(locate, patterns, library, cutdirs=0):
	files = []
	path = ""
	for pattern in patterns:
		files = files + os.popen(locate + " " + pattern).readlines()
	if len(files) == 0:
		print "ERROR:  " + library + " not found!"
	if len(files) == 1:
		path = pathTruncated(files[0],cutdirs)
		print library + " found:\t" + path
	if len(files) > 1:
		paths = []
		for filename in files:
			paths.append(pathTruncated(filename,cutdirs))
		uniquepaths = unique(paths)
		if len(uniquepaths) >= 1:
			path = pathTruncated(files[0],cutdirs)
		if len(uniquepaths) > 1:
			print "WARNING: " + library + " directory ambiguous:"
			for uniquepath in uniquepaths:
				print "\t " + uniquepath
			print "\t Using " + path + " as directory for " + library
		print library + " found:\t" + path
	return path
	
# search directories for headers and libraries needed for ShapeLens++ compilation
# stores the results in filename
def createConfig(filename):
	print "@@@ Searching for include and library directories @@@"
	# set up necessary paths (and searchkeys)
	paths = {}
	paths["SHAPELENSPATH"] = ""
	paths["NUMLAPATH"] = ""
	paths["GSLINCLPATH"] = ""
	paths["GSLLIBPATH"] = ""
	paths["BOOSTINCLPATH"] = ""
	paths["BOOSTSANDBOXINCLPATH"] = ""
	paths["ATLASINCLPATH"] = ""
	paths["ATLASLIBPATH"] = ""
	paths["CBLASLIBPATH"] = ""
	paths["G2CLIBPATH"] = ""
	paths["CFITSIOINCLPATH"] = ""
	paths["CFITSIOLIBPATH"] = ""
	# check for locate on the system
	locate = os.popen("which locate").readlines()
	if len(locate) == 1:
		locate = locate[0].strip()
		paths["SHAPELENSPATH"] = checkWithLocate(locate, ["shapelens/include/ShapeLens.h"], "ShapeLens++ headers",1)
		paths["NUMLAPATH"] = checkWithLocate(locate, ["numla/NumMatrix.h"], "numla headers")
		paths["GSLINCLPATH"] = checkWithLocate(locate, ["gsl/gsl_math.h"], "GSL headers",1)
		paths["GSLLIBPATH"] = checkWithLocate(locate, ["libgsl.so","libgsl.a"], "GSL library")
		paths["BOOSTINCLPATH"] = checkWithLocate(locate, ["boost/config.hpp"],"Boost headers",1)
		paths["BOOSTSANDBOXINCLPATH"] = checkWithLocate(locate, ["boost/numeric/bindings/traits/ublas_matrix.hpp"],"Boost sandbox headers",4)
		paths["ATLASINCLPATH"] = checkWithLocate(locate, ["atlas_lapack.h"],"ATLAS headers")
		paths["ATLASLIBPATH"] = checkWithLocate(locate, ["libatlas.a","libatlas.so"],"ATLAS library")
		paths["CBLASLIBPATH"] = checkWithLocate(locate, ["libcblas.a","libcblasb.so"],"CBLAS library")
		paths["G2CLIBPATH"] = checkWithLocate(locate, ["libg2c.a","libg2c.so"],"G2C library")
		paths["CFITSIOINCLPATH"] = checkWithLocate(locate, ["fitsio.h"],"cfitsio headers")
		paths["CFITSIOLIBPATH"] = checkWithLocate(locate, ["libcfitsio.a","libcfitsio.so"],"cfitsio library")
		print "\n\nPlease check the directories for headers and libraries in the file " + filename + "\n"
	else:
		print "\n\nLocate does not exist on this system."
		print "Please open " + filename + " and specify the directories for headers and libraries." + "\n"
	configfile = open(filename,"w")
	keys = paths.keys()
	keys.sort()
	for key in keys:
		configfile.write(key+"\t"+paths[key]+"\n")
	configfile.close()
	sys.exit(1)

# open file created by searchPaths()
# does some additional error checking
def openConfig(filename, numentries):
	print "@@@ Opening configuration file " + filename + " @@@"
	configfile = open(filename)
	lines = configfile.readlines()
	configfile.close()
	paths = {}
	for line in lines:
		# exclude potential comment lines
		if line.find("#",0,0) == -1:
			chunks = line.split("\t")
			if len(chunks) == 2:
				key = chunks[0]
				value = chunks[1].strip()
				# check if directory exists
				# because of $(HOME) replacement in pathTruncated(), we have to re-replace
				if os.path.isdir(value.replace("$(HOME)",os.environ["HOME"])):
					paths[key] = value
				else:
					print "ERROR:  " + chunks[0] + " directory " + chunks[1].strip() + " does not exist"
					print "\tPlease correct this entry in the config file " + filename
					sys.exit(1)
			else:
				print "ERROR:   config file format is: <NAME><TAB><DIRECTORY>"
				print "\tPlease correct the config file " + filename
				sys.exit(1)
	if len(paths) != numentries:
		print "ERROR:   Number of entries in config file does not match requirement of " + str(numentries)
		sys.exit(1)
	return paths


# first of all: create/open the configuration file, which lists the include and lib directories
if os.path.isfile(CONFIGFILE):
	config = openConfig(CONFIGFILE,12)
else:
	createConfig(CONFIGFILE)

# insert PATH and CFLAGS in Makefile
print "@@@ Creating new Makefile @@@"
infile = open('Makefile.in','r') 
makefile = infile.read()
infile.close()
paths = "INSTALLPATH = "+INSTALLPATH+"\nSHAPELENSPATH = "+config["SHAPELENSPATH"]+"\nNUMLAPATH = "+config["NUMLAPATH"]

targets = "all: lib"
doctarget = ""
if (SHARED == True):
        targets = targets + " shared"
if (DOCUMENTATION == True):
	targets = targets + " docs"
	doctarget = "docs: $(HEADERS)\n\tdoxygen Doxyfile"
targets = targets + " progs"

# loop thru paths and append avery occurence of "LIBPATH" to otherlib (same for "INCLPATH")
otherincl = []
otherlib = []
for key, value in config.iteritems():
	if key.find("INCLPATH") != -1:
		otherincl.append(value)
	if key.find("LIBPATH") != -1:
		otherlib.append(value)
otherincl = unique(otherincl)
otherlib = unique(otherlib)
if len(otherincl) >=1:
	otherincl = " -I" + " -I".join(otherincl)
if len(otherlib) >=1:
	otherlib = " -L" + " -L".join(otherlib)

makefile = makefile.replace("???PATHS???",paths)
makefile = makefile.replace("???SUBDIR???",SUBDIR)
makefile = makefile.replace("???COMPILER???",COMPILER)
makefile = makefile.replace("???CFLAGS???",FLAGS)
makefile = makefile.replace("???OTHERINCL???",otherincl)
makefile = makefile.replace("???OTHERLIB???",otherlib)
makefile = makefile.replace("???LIBS???",LIBS)
makefile = makefile.replace("???TARGETS???",targets)
makefile = makefile.replace("???DOCTARGET???",doctarget)	

outfile = open('Makefile','w')
outfile.write(makefile)
outfile.close()

# insert INPUT in Doxyfile
if (DOCUMENTATION == True):
	print "@@@ Creating new Doxyfile @@@"
	infile = open('Doxyfile.in','r')
	doxyfile = infile.read()
	infile.close()

	input = config["NUMLAPATH"]+" "+config["SHAPELENSPATH"]+"/include"
	doxyfile = doxyfile.replace("???INPUT???",input)

	outfile = open('Doxyfile','w')
	outfile.write(doxyfile)
	outfile.close()

# create directories
if not os.path.isdir("progs/"+SUBDIR):
	print "@@@ Creating directory progs/"+SUBDIR + " @@@"
	os.system("mkdir -p progs/"+SUBDIR)
if not os.path.isdir("lib/"+SUBDIR):
	print "@@@ Creating directory lib/"+SUBDIR + " @@@"
	os.system("mkdir -p lib/"+SUBDIR)
if DOCUMENTATION == True and not os.path.isdir("docs"):
	print "@@@ Creating directory docs @@@"
	os.system("mkdir -p docs")

# check environment variable PATH and LD_LIBRARY_PATH
pwd = os.getcwd()
if os.environ["PATH"].find(pwd +"/progs/" + SUBDIR) == -1:
	print "CAUTION: "+ pwd +"/progs/" + SUBDIR + " not in $PATH"
if SHARED== True and os.environ["LD_LIBRARY_PATH"].find(pwd +"/lib/" + SUBDIR) == -1:
	print "CAUTION: "+ pwd +"/lib/" + SUBDIR + " not in $LD_LIBRARY_PATH"
	print "\t Please set the environment variable LD_LIBRARY_PATH for using the shared library"
