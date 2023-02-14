# CLBLM Makefile v 0.9
# Author: Chris Brodowski (cbrodows@aer.com)
#
# This Makefile has been restructured so that CLBLM is compiled as a series of
# libraries.  These then get linked in when the final executable is built.
# The libraries are:
# -libCommon (common routines used by CLBLM)
# -libClbl (CLBLM src)
# -libScene (Scene tool)
# -libJson (JSON parser library)

CLBLM_DIR = src/clblm_src
SCENE_DIR = src/scene_tool
COMMON_DIR = src/common
APP_DIR = src/app
JSON_DIR = src/json
MODULES_DIR = src/modules/

# ifort or gfortran; *** CON NOT HAVE SPACE BEFORE OR AFTER ***
CL?=ifort

LIB_TYPE = "static"	# static or dynamic library linking; this is a TODO

# GFORTRAN setup
# NOTE: As of right now, CLBLM will not build with gfortran using -std=f2003 or -std=f2008, since not all the source code actually adheres to either of these standards.
ifeq ($(CL),gfortran)
    #CFLAGS = -fbounds-check -fdefault-real-8 -frecord-marker=4 -g -fpic -DGFORTRAN -cpp -ffree-line-length-none -fno-align-commons -J ./modules/
    CFLAGS = -fno-range-check -fbounds-check -fdefault-real-8 -frecord-marker=4 -g -fpic -DGFORTRAN -cpp -ffree-line-length-none -fno-align-commons -J $(MODULES_DIR)
endif

# Intel Compiler setup
ifeq ($(CL),ifort)
    CFLAGS = -i8 -r8 -module $(MODULES_DIR)
endif

# TODO: Portland Group setup


NETCDF_LIBRARY_LINK = $(shell nc-config --flibs)
NETCDF_INCLUDES = $(shell nc-config --fflags)
#NETCDF_INCLUDES = $(shell nc-config --fflags) -I/usr/lib64/gfortran/modules


# Common Library
COMMON_SRC = $(COMMON_DIR)/util_linux_intel.f90 \
	          $(COMMON_DIR)/planet_earth.f90 \
	          $(COMMON_DIR)/clblm_Util.f90 \
	          $(COMMON_DIR)/clblm_ConstParam.f90

COMMON_OBJ = $(patsubst %.f90,%.o,$(COMMON_SRC))
COMMON_LIB: $(COMMON_OBJ)
	ar rv libCommon.a $(COMMON_OBJ)


# "src"/CLBLM library
CLBLM_SRC = $(CLBLM_DIR)/clblm_Config.f90 \
	         $(CLBLM_DIR)/clblm_Spectrum.f90 \
	         $(CLBLM_DIR)/clblm_EMLAY.f90 \
	         $(CLBLM_DIR)/clblm_FileIO.f90 \
	         $(CLBLM_DIR)/clblm_Scene.f90 \
	         $(CLBLM_DIR)/clblm_AtmPath.f90 \
	         $(CLBLM_DIR)/clblm_LineData.f90 \
	         $(CLBLM_DIR)/clblm_setDV.f90 \
	         $(CLBLM_DIR)/clblm_ScanFilter.f90 \
	         $(CLBLM_DIR)/clblm_FFT.f90 \
	         $(CLBLM_DIR)/clblm_PostProc.f90 \
	         $(CLBLM_DIR)/clblm_ODLAY_TIPS.f90 \
	         $(CLBLM_DIR)/clblm_ODLAY_CommonSub.f90 \
	         $(CLBLM_DIR)/clblm_ODLAY_XSect.f90 \
			 $(CLBLM_DIR)/read_module.f90 \
			 $(CLBLM_DIR)/mt_ckd_h2o_module.f90 \
	         $(CLBLM_DIR)/clblm_ODLAY_Continuum.f90 \
	         $(CLBLM_DIR)/clblm_ODLAY_LineF4.f90 \
	         $(CLBLM_DIR)/clblm_ODLAY.f90 \
            $(CLBLM_DIR)/solar_cycle.f90\
            $(CLBLM_DIR)/clblm_Solar.f90\
	         $(CLBLM_DIR)/clblm_noScattRT.f90 \
	         $(CLBLM_DIR)/clblm_noScattJac.f90 \
	         $(CLBLM_DIR)/clblm_Drivers.f90 \
	         $(CLBLM_DIR)/clblm_Modules.f90

CLBLM_OBJ = $(patsubst %.f90,%.o,$(CLBLM_SRC))
CLBLM_LIB: $(CLBLM_OBJ)
	ar rv libClbl.a $(CLBLM_OBJ)


# "src"/Scene tool library
SCENE_SRC = $(SCENE_DIR)/scene_Config.f90\
            $(SCENE_DIR)/scene_UnitConv.f90\
            $(SCENE_DIR)/scene_TAPE5.f90\
            $(SCENE_DIR)/scene_ModelAtm.f90\
            $(SCENE_DIR)/scene_Surface.f90\
            $(SCENE_DIR)/scene_Profile.f90\
            $(SCENE_DIR)/scene_Geometry.f90\
            $(SCENE_DIR)/scene_netcdfWriter.f90\
            $(SCENE_DIR)/scene_Modules.f90

SCENE_OBJ = $(patsubst %.f90,%.o,$(SCENE_SRC))
SCENE_LIB: $(SCENE_OBJ)
	ar rv libScene.a $(SCENE_OBJ)


# JSON library
JSON_SRC = $(JSON_DIR)/jsonlib/jsonModule.f90 \
	        $(JSON_DIR)/jsonlib/jsonFile.f90 \
	        $(JSON_DIR)/DataTypes.f90 \
	        $(JSON_DIR)/validateJSONmodule.f90 \
	        $(JSON_DIR)/StringUtilities.f90 \
	        $(JSON_DIR)/JsonConfig.f90

JSON_OBJ = $(patsubst %.f90,%.o,$(JSON_SRC))
JSON_LIB: $(JSON_OBJ)
	ar rv libJson.a $(JSON_OBJ)


# Can do a list of executables here, as long as you watch what you link, you
# can compile all of them in one fell swoop
EXEC_SRC = $(APP_DIR)/CheckCompile.f90 \
		     $(APP_DIR)/scene_writer.f90 \
		     $(APP_DIR)/build_solar.f90 \
		     $(APP_DIR)/clblm_main.f90

EXEC_OBJ = $(patsubst %.f90,%.o,$(EXEC_SRC))


# ----------------------- Build/Link Rules -----------------------

# Include all .mod and NetCDF stuff
INCLUDES = -I$(MODULES_DIR) $(NETCDF_INCLUDES)

# NOTE: Repeating the linking command a second time eliminates cross-dependency
# 	issues, and the compiler will ignore the duplicate/extra links
# 	once successfully ordered.
LINKING_COMMAND = -L./ -lJson -lCommon -lClbl -lScene -lJson -lCommon -lClbl -lScene $(NETCDF_LIBRARY_LINK)

# This builds all libraries.
libs: JSON_LIB COMMON_LIB SCENE_LIB CLBLM_LIB

# This is just to see if the linking works with JSON and NETCDF
compileTest: libs $(EXEC_OBJ)
	$(CL) $(CFLAGS) $(INCLUDES) -o compileTest $(APP_DIR)/CheckCompile.o $(LINKING_COMMAND)

# scene writer
scene_writer: libs $(EXEC_OBJ)
	$(CL) $(CFLAGS) $(INCLUDES) -o scene_writer $(APP_DIR)/scene_writer.o $(LINKING_COMMAND)

# NetCDF file builder
build_solar: libs $(EXEC_OBJ)
	$(CL) $(CFLAGS) $(INCLUDES) -o build_solar $(APP_DIR)/build_solar.o $(LINKING_COMMAND)

# clblm main
clblm: libs $(EXEC_OBJ)
	$(CL) $(CFLAGS) $(INCLUDES) -o clblm $(APP_DIR)/clblm_main.o $(LINKING_COMMAND)

#
all: clblm scene_writer build_solar

# ------------------------ Compile Rules ------------------------------

.f90.o:
	$(CL) -c $(CFLAGS) $(INCLUDES) $< -o $@

.PHONY: clean

clean-compile:
	@find . -iname '*.o' -delete
	@find . -iname '*.mod' -delete

# TODO: add clean of executables to this directory
clean:
	$(MAKE) clean-compile
	@rm -f *.a

.SUFFIXES: .o .f90


