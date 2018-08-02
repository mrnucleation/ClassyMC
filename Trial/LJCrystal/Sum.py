infile = open("fort.2", "r")

suma = 0.0
for line in infile:
    col = line.split()
    suma += float(col[-1])

print(suma)

