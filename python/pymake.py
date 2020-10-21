#!/usr/local/bin/python

import distutils.sysconfig
import string, sys

configopts = {}

maketemplate = """
PYLIB=%(pythonlib)s
PYINC=-I%(pythoninc)s
LIBS=%(pylibs)s
OPTS=%(pyopt)s
SRC_MAIN := $(PYTHON)/forpy_mod.f90 $(PYTHON)/Analysis_Python.f90  $(PYTHON)/Python_CommonTypes.f90 $(PYTHON)/Sim_Python.f90 $(SRC_MAIN)
PACKAGE_FLAGS += -DEMBPYTHON   $(PYINC) $(LIBS)
OBJ_LIBRARY += $(PYLIB)  


$(OBJ)/Python_CommonTypes.o:  $(OBJ)/Common_BoxData.o $(OBJ)/Box_CubicBox.o $(OBJ)/Box_OrthoBox.o
$(OBJ)/Sim_Python.o:  $(OBJ)/Sim_MonteCarlo.o
$(OBJ)/Common.o += $(OBJ)/forpy_mod.o $(OBJ)/Sim_Python.o

$(OBJ)/Analysis_Python.o: $(OBJ)/forpy_mod.o  $(OBJ)/Template_Analysis.o $(OBJ)/Input_Format.o $(OBJ)/Common_Analysis.o  $(OBJ)/Common_BoxData.o  $(OBJ)/Python_CommonTypes.o
"""

configopts['pythonlib'] = distutils.sysconfig.get_config_var('LIBPL') \
        + '/' + \
        distutils.sysconfig.get_config_var('LIBRARY')
configopts['pythoninc'] = ''
configopts['pylibs'] = ''

#for dir in string.split(distutils.sysconfig.get_config_var('INCLDIRSTOMAKE')):
for dir in distutils.sysconfig.get_config_var('INCLDIRSTOMAKE').split():
    configopts['pythoninc'] += '-I%s ' % (dir,)
#for dir in string.split(distutils.sysconfig.get_config_var('LIBDIR')):
for dir in distutils.sysconfig.get_config_var('LIBDIR').split():
    configopts['pylibs'] += '-L%s ' % (dir,)

configopts['pylibs'] += distutils.sysconfig.get_config_var('MODLIBS') \
        + ' ' + \
        distutils.sysconfig.get_config_var('LIBS') \
        + ' ' + \
        distutils.sysconfig.get_config_var('SYSLIBS')

configopts['pyopt'] = distutils.sysconfig.get_config_var('OPT')

targets = ''

for arg in sys.argv[1:]:
    targets += arg + ' '

configopts['programs'] = targets
with open("../Python.Makefile", "w") as outfile:
    outfile.write("{}\n".format(maketemplate % configopts))
#    for arg in sys.argv[1:]:
#        outfile.write("{}\n".format("%s: %s.o\n\tgcc %s.o $(LIBS) $(PYLIB) -o %s" % (arg, arg, arg, arg))
#        outfile.write("{}\n".format("%s.o: %s.c\n\tgcc %s.c -c $(PYINC) $(OPTS)" % (arg, arg, arg))
#
#    print( "clean:\n\trm -f $(PROGRAMS) *.o *.pyc core")
