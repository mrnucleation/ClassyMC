
#==================================================
'''
 Script designed to convert XSF format over to Classy's Data file.
'''

#==================================================
class BoxTypeError(Exception):
    pass

#==================================================
def xsftoclassy(inname, outname):
    with open(inname, "r") as infile:
        header, coords = GetXSFFrame(infile)
    for line in coords:
        print(line)
    

    boxtype, cellvectors = GetBoxType(header)

    nAtoms = {}
    typelist = {}
    typecounter = 0
    for line in coords:
        atom = line.split()
        atmtype = atom[0]
        if atmtype not in typelist:
            typecounter += 1
            typelist[atmtype] = typecounter
            nAtoms[atmtype] = 1
        else:
            nAtoms[atmtype] += 1


    with open(outname, "w") as outfile:
        outfile.write("boxtype %s\n"%(boxtype))
        if boxtype is "nobox":
            aoutfile.write("\n")
        elif boxtype is "ortho":
            outvec = [str(cellvectors[i][i]) for i in range(3)]
            outfile.write("dimension %s %s %s\n"%(tuple(outvec)))
        elif boxtype is "triclinic":
            raise BoxTypeError("Trilcinic has not been implimented in Classy yet!")
        else:
            raise BoxTypeError("Unknown Box Type Found")
        outstr = ' '.join([str(nAtoms[key]) for key in nAtoms])
        outfile.write("molmin   %s\n"%(outstr))
        outfile.write("molmax   %s\n"%(outstr))
        outfile.write("mol   %s\n"%(outstr))
        label = {}
        for i in range(1,typecounter+1):
            label[i] = 0
        for line in coords:
            atom = line.split()
            atmtype, x, y, z = tuple(atom)
            outtype = typelist[atmtype]
            label[outtype] += 1
            outfile.write("%s  %s  %s  %s  %s  %s\n"%(outtype, label[outtype], 1, x, y, z))

    print("Finished Converting %s to %s"%(inname, outname))





#==================================================
def GetBoxType(header):
    from itertools import combinations
    from math import fabs

    #Check to see if the coodinates even require a simulation box or are given
    #as a cluster configuration by searching the header for the proper keywords
    for i, line in enumerate(header):
        if "PRIMVEC" in line:
            startline = i
            break
             
        elif "ATOMS" in line:
            return "nobox", None

    #If a box is found, collect the cell vectors
    vectors = [ [0.0 for x in range(3)] for y in range(3)]
    for i, line in enumerate(header[startline+1:startline+1+3]):
        col = line.split()
        vectors[i][0] = float(col[0])
        vectors[i][1] = float(col[1])
        vectors[i][2] = float(col[2])

    #Check of the off diagonal vectors are non-zero which indicates a triclinic box
    for i,j in combinations(list(range(3)), 2):
        if i == j:
            continue
        if vectors[i][j] > 1e-6:
            return "triclinic", vectors
        
    #If the vectors aren't triclinic they should be orthogonal. 
    #If it isn't......there's something rotten in the state of Denmark. 
    return "ortho", vectors       
#==================================================
def GetXSFFrame(infile):
    last_pos = infile.tell()
    line = infile.readline()
#    if line == '':
#        raise IOError
    col = line.split()

    header = []
    header.append(line)


    coords = []
 
    blankCnt = 0
    while True: 
        last_pos = infile.tell()
        line = infile.readline().strip()
        if line == '' or (line is None) or (len(line.split()) < 1):
            blankCnt += 1
            if blankCnt > 3:
                return header, coords
        else:
            blankCnt = 0
#        print line
        col = line.split()
        if len(col) < 4:
            header.append(line)
        else:
            try:
                for item in col[1:4]:
                    float(item)
            except TypeError:
                header.append(line)

            finally:
                coords.append(line)


        
    return header, coords


if __name__ == "__main__":
    import sys
    inname = sys.argv[1]
    outname = sys.argv[2]
    xsftoclassy(inname, outname)



