#===================================================
#Python Bindings for the ClassyMC code base. 
#===================================================
import os
from ctypes import c_int
import cffi

#===================================================
class ClassyCrash(Exception):
    pass
#===================================================
class ClassyMC(object):
    #------------------------------------
    def __init__(self, startscript=None, verbose=False):
        self.FFI = cffi.FFI()
        self.ClassyLib = self.__startlibrary()

        self.verbose = verbose
        self.prologue = False

        self.movecount = 0

        # Classy Library Function Signatures
        self.FFI.cdef("""void Classy_ReadScript(int strlen, char *str);""")
        self.FFI.cdef("""void Classy_RunMove();""")
        self.FFI.cdef("""void Classy_FullSim();""")
        self.FFI.cdef("""void Classy_SetLogfile();""")
        self.FFI.cdef("""void Classy_RunPrologue();""")
        self.FFI.cdef("""void Classy_ScreenOut(long cyclenum);""")
        self.FFI.cdef("""void Classy_ForceMove(int movenum);""")

        if startscript is not None:
            self.readscript(startscript)

    #------------------------------------
    def __startlibrary(self):
        return self.FFI.dlopen("./libclassymc.so")

    #------------------------------------
    def readscript(self, infile):
        '''
        # Passes the name of an input file into Classy's Script engine
        # in the same fashion as if the script was passed in via the command
        # line. Can be used to initialize Classy in preparation for Python
        # level control.
        '''
        if self.verbose:
            print("Classy: Reading from file %s"%(infile))
        c_infile = self.FFI.new("char[]", infile)
        c_insize = c_int()
        c_insize = len(infile)
        self.ClassyLib.Classy_ReadScript(c_insize, c_infile)
    #------------------------------------
    def setlogfile(self, filename):
        '''
        # Passes the name of an input file into Classy's Script engine
        # in the same fashion as if the script was passed in via the command
        # line. Can be used to initialize Classy in preparation for Python
        # level control.
        '''
        if self.verbose:
            print("Classy: Setting Log File to file %s"%(filename))
        c_infile = self.FFI.new("char[]", filename)
        c_insize = c_int()
        c_insize = len(infile)
        self.ClassyLib.Classy_SetLogFile(c_insize, c_infile)

    #------------------------------------
    def runsimulation(self):
        '''
        # Tells Classy to perform an entire simulation. Equivalent
        # to the "run" command in a standard Classy script
        '''
        self.ClassyLib.Classy_FullSim()
        self.prologue = True

    #------------------------------------
    def runmove(self):
        '''
        # Tells Classy to perform a single randomly selected
        # MC move and returns back to Python
        '''
        if not self.prologue:
            self.ClassyLib.Classy_RunPrologue()
            self.prologue = True
        self.ClassyLib.Classy_RunMove()
        self.movecount += 1

    #------------------------------------
    def forcemove(self, movenum):
        pass
        '''
        # Tells Classy to perform a specified Monte Carlo move. Useful
        # if the move selection is handled at the Python level.
        '''
        self.ClassyLib.Classy_ForceMove()
        self.movecount += 1
    #------------------------------------
    def printstat(self):
        '''
        # Tells Classy to print mid-simulation screen output
        '''
        self.ClassyLib.Classy_ScreenOut(self.movecount)
    #------------------------------------
    def getpositions(self, boxnum):
        '''
        # Collects the atomic positions of a given simulation box from
        # 
        '''
        c_boxnum = c_int()
        c_boxnum = boxnum
        pass
    #------------------------------------
    def feedpositions(self):
        pass
    #------------------------------------
    def feedenergy(self, E_New):
        pass

    #------------------------------------


#===================================================
