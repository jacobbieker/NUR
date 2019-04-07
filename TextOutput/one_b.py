import numpy as np
import matplotlib.pyplot as plt

def random_generator(seed, m=2 ** 64 - 1, a=2349543, c=913842, a1=21, a2=35, a3=4, a4=4294957665):
    """
        Generates psuedorandom numbers with a combination of (M)LCC, 64 bit shift, and MWC
    :param seed: Seed to use
    :param m: Determines period of the MLCC
    :param a: For the MLCC
    :param c: For the MLCC
    :param a1: For the first bit shift
    :param a2: For the second bit shift
    :param a3: For the third bit shift
    :param a4: For the MWC
    :return:
    """

    # First linear congruential generator
    # While true, so the generator never stops making new numbers
    # This is used to make sure teh XOR shift is 64 bit
    bit_64 = 0xffffffffffffffff
    while True:
        # This is MLCC part
        generated_number = (a * seed + c) % m
        # Now bit shift
        generated_number = generated_number ^ (generated_number >> a1) & bit_64
        generated_number = generated_number ^ (generated_number << a2) & bit_64
        generated_number = generated_number ^ (generated_number >> a3) & bit_64

        # Now MWC part
        mwc_out = a4 * (generated_number & (2 ** 32 - 1)) + (generated_number >> 32)

        seed = mwc_out
        # m = third_output
        mwc_out = mwc_out / m

        if mwc_out > 1.:
            # Have to make it between 1 and 0, so mod 1. makes sure its between 0 and 1 now
            close_to_final = mwc_out % 1.
        else:
            close_to_final = mwc_out

        yield close_to_final

def one_b(rand_gen):
    first_thousand = []
    first_thousand_x_1 = [0.0]
    # 1b first one
    for i in range(1000):
        # This is x_i+1
        first_thousand.append(next(rand_gen))
        if i > 0:
            # This is x_i
            first_thousand_x_1.append(first_thousand[i - 1])

    # Now plot xi+1 vs xi

    plt.scatter(first_thousand, first_thousand_x_1)
    plt.xlabel("$X_i$")
    plt.ylabel("$X_{i+1}$")
    plt.title("$X_{i+1}$ vs. $X_i$")
    plt.savefig("./plots/Xi_Xi_1.png", dpi=300)
    plt.cla()

    first_million = []
    for i in range(1000000):
        first_million.append(next(rand_gen))

    plt.hist(first_million, bins=np.linspace(0.0, 1.0, 20))
    plt.savefig("./plots/1000000_rand.png", dpi=300)
    plt.cla()
