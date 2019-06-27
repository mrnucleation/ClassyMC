#===================================================
#Python Bindings for the ClassyMC code base. 
#===================================================
import os
from ctypes import c_int
import cffi


#===================================================
class ClassyMC(object):
    #------------------------------------
    def __init__(self, startscript=None):
        self.FFI = cffi.FFI()
        self.ClassyLib = self.__startlibrary()

        # Classy Library Function Signatures
        self.FFI.cdef("""void Classy_ReadScript(int strlen, char *str);""")
        self.FFI.cdef("""void Classy_RunMove();""")
        self.FFI.cdef("""void Classy_FullSim();""")
#        self.FFI.cdef("""void Classy_RunMove();""")

        if startscript is None:
            self.readscript(startscript)

    #------------------------------------
    def __startlibrary(self):
        return FFI.dlopen("./libclassymc.so")

    #------------------------------------
    def readscript(self, infile):
        '''
        # Passes the name of an input file into Classy's Script engine
        # in the same fashion as if the script was passed in via the command
        # line.
        '''
        c_infile = FFI.new("char[]", infile)
        c_insize = c_int()
        c_insize = len(infile)
        self.ClassyLib.Classy_ReadScript(c_insize, c_infile)
    #------------------------------------
    def runsimulation(self):
        pass
        '''
        # Tells Classy to perform an entire simulation. Equivalent
        # to the "run" command in a standard Classy script
        '''
        self.ClassyLib.Classy_FullSim()

    #------------------------------------
    def runmove(self):
        '''
        # Tells Classy to perform a single randomly selected
        # MC move and returns back to Python
        '''
        self.ClassyLib.Classy_RunMove()
    #------------------------------------
    def forcemove(self, movenum):
        pass
        '''
        # Tells Classy to perform a specified Monte Carlo move.
        '''
#        self.ClassyLib.Classy_ReadScript()
    #------------------------------------
    def getpositions(self):
        pass
    #------------------------------------


#===================================================
