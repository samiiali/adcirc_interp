################################################################################
### Makefile template for user ESMF application, leveraging esmf.mk mechanism ##
################################################################################

################################################################################
### Finding and including esmf.mk ##############################################

# Note: This fully portable Makefile template depends on finding environment
#       variable "ESMFMKFILE" set to point to the appropriate "esmf.mk" file,
#       as is discussed in the User's Guide.
#       However, you can still use this Makefile template even if the person
#       that installed ESMF on your system did not provide for a mechanism to
#       automatically set the environment variable "ESMFMKFILE". In this case
#       either manually set "ESMFMKFILE" in your environment or hard code the
#       location of "esmf.mk" into the include statement below.
#       Notice that the latter approach has negative impact on portability.

# Project name and version
Proj := Test
Version := Debug

#paths for Project (Ppath) Object files (Opath) and binary path (Bpath)
Ppath := /workspace/ADCIRC_Interpolation/test2/ali_test/${Proj}
Opath := .
Bpath := .

ESMFMKFILE=/workspace/Library/ESMF/install/lib/libg/Linux.gfortran.64.mpich2.default/esmf.mk

include $(ESMFMKFILE)

COMPILE = mpif90


Debug: all
Release: all

all: $(Proj)


################################################################################
### Compiler and linker rules using ESMF_variables supplied by esmf.mk ########


%.o: %.f90
	
	$(COMPILE) -g -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) \
          $(ESMF_F90COMPILEFREENOCPP) $<

%.o: %.F90
	$(COMPILE) -g -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) \
          $(ESMF_F90COMPILEFREENOCPP) $<

%: %.o
	$(COMPILE) -g $(ESMF_F90LINKOPTS) $(ESMF_F90LINKPATHS) \
          $(ESMF_F90LINKRPATHS) -o $@ $^ $(ESMF_F90ESMFLINKLIBS)
