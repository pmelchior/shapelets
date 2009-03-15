import os

# check whether CODEDIR and ITALIBSINCLPATH are define
CODEDIR = ""
ITALIBSPATH = ""
ITALIBSINCLPATH = ""
ITALIBSLIBPATH = ""
if os.environ.get("CODEDIR") == None:
	print "CODEDIR variable not defined"
	os.exit(1)
else:
	CODEDIR = os.environ.get("CODEDIR")
if os.environ.get("ITALIBSPATH") == None:
	print "ITALIBSPATH variable not defined"
	os.exit(1)
else:
	ITALIBSPATH = os.environ.get("ITALIBSPATH")
	ITALIBSINCLPATH = ITALIBSPATH + "/include"
	ITALIBSLIBPATH = ITALIBSPATH + "/lib"
	ITALIBSPROGPATH = ITALIBSPATH + "/bin"
	if os.environ.get("SUBDIR"):
		ITALIBSLIBPATH = ITALIBSLIBPATH + "/" + os.environ.get("SUBDIR")
		ITALIBSPROGPATH = ITALIBSPROGPATH + "/" + os.environ.get("SUBDIR")
	if os.path.exists(ITALIBSINCLPATH) == False:
		os.system("mkdir " + ITALIBSINCLPATH)
	if os.path.exists(ITALIBSLIBPATH) == False:
		os.system("mkdir " + ITALIBSLIBPATH)
	if os.path.exists(ITALIBSPROGPATH) == False:
		os.system("mkdir " + ITALIBSPROGPATH)

# check for numla files
if os.path.isdir(CODEDIR+"/numla"):
	if os.path.exists(ITALIBSINCLPATH + "/numla") == False:
		os.system("ln -s " + CODEDIR+"/numla " + ITALIBSINCLPATH + "/numla");
		print "Check that libatlas.so, liblapack_atlas.so, libcblas.so and liblapack.so are either in standard locations or in " + ITALIBSLIBPATH
else:
	print CODEDIR+"/numla does not exist."
	os.exit(1)

# add shapelens files
if os.path.isdir(CODEDIR + "/shapelens"):
	if os.path.exists(ITALIBSINCLPATH + "/shapelens") == False:
		os.system("ln -s " + CODEDIR+"/shapelens/include " + ITALIBSINCLPATH + "/shapelens")
	if os.path.exists(ITALIBSLIBPATH+"/libshapelens.so") == False:
		os.system("ln -s " + CODEDIR+"/shapelens/lib/libshapelens.so "+ITALIBSLIBPATH+"/")
else:
	print CODEDIR + "/shapelens does not exist."
	os.exit(1)
