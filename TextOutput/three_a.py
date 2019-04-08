import numpy as np

"""

Dont write down formula

Produce set of random numbers

Returning log likliehood corresponding to all those numbers

Physcially write down log-likliehood for satallite equation and unknown Nsat

Disregard all terms without dependence on a, b,c, or A

Then do a,b,c that maxes it


use 1000 numbers from 2 to use as y in the formula

"""

# Read in file, one at a time for memory limits

def three_a():
    files = ["satgals_m11.txt", "satgals_m12.txt", "satgals_m13.txt", "satgals_m14.txt", "satgals_m15.txt"]

    for filename in files:
        with open(filename, "r") as mass_haloes:
            mass_haloes.readline()
            mass_haloes.readline()
            mass_haloes.readline()
            # Now next one will be number of haloes
            num_haloes = mass_haloes.readline()
            haloes = []
            for line in mass_haloes:
                if "#" in line:
                    # new halo
                    haloes.append([])
                else:
                    try:
                        radius = float(line.split(" ")[0])
                        haloes[-1].append(radius)
                    except:
                        print("No Radius, but not #")
                        print(line)
