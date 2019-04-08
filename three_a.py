import numpy as np
import matplotlib.pyplot as plt

"""

From before, changing the random seed changed the distribution of satellites with regards to x

So the log likelihood we need is the a,b,c that puts the peak of the distribution over the peak of the 
distribution in each halo mass bin text file



"""

# Read in file, one at a time for memory limits

def three_a():

    files = ["satgals_m11.txt", "satgals_m12.txt", "satgals_m13.txt", "satgals_m14.txt", "satgals_m15.txt"]

    for filename in list(reversed(files)):
        summed_radii = []
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
                        summed_radii.append(radius)
                    except:
                        print("No Radius, but not #")
                        print(line)

            haloes = np.asarray(haloes)
            summed_radii = np.asarray(summed_radii)
            bins = np.linspace(0,2.5,101)
            plt.hist(summed_radii, bins=bins)
            plt.yscale("log")
            plt.xlabel("Radius in X")
            plt.title(filename)
            plt.savefig("./plots/{}.png".format(filename.split(".")[0]))
            plt.cla()
