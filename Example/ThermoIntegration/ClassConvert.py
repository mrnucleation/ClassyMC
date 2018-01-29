
from math import sqrt
infile = open("FCC.xyz", "r")
outfile = open("FCC.clssy", "w")

cnt = 0
scale = 2.0**(2.0/3.0)
outfile.write( "boxtype cube \n" )
outfile.write( " dimension " +  str(8.0*scale) + "\n" )

for line in infile:
    try:
        col = line.split()
        x = (float(col[1]) - 4.0) * scale
        y = (float(col[2]) - 4.0) * scale
        z = (float(col[3]) - 4.0) * scale

        cnt += 1
        outfile.write(' '.join([str(x) for x in [1, cnt, 1, x, y, z, "\n"]]) )
    except:
        continue
