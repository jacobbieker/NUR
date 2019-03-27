import numpy as np
import matplotlib.pyplot as plt


def poisson(lam, k):
    # TODO Makes almost zero for some values, need to change calculation
    """
    Possion Function (lambda^k)*(e^-lambda)/(k!)
    :param lam:
    :param k:
    :return:
    """
    exponent_part = np.float64(np.exp(-1. * lam))
    first_part = np.float64(lam ** k)
    # Factorial part....
    fact = np.float64(1)
    for i in range(1, np.int64(k + 1)):
        fact *= i

    final_result = np.float64((first_part * exponent_part) / fact)
    return final_result


print(poisson(1, 0))
print(poisson(5, 10))
print(poisson(3, 21))
print(poisson(2.6, 40))

seed = 1337

"""

Add x^2*4*pi for a 1D integral to make it the spherical one

Use midpoint rule, since there is a 0

A = 1/everything else


"""


def random_generator(seed, m=2**64, a=2349543, c=913842, a1=21, a2=35, a3=4, a4=4294957665):
    """
    Generates psuedorandom numbers with a combination of (M)LCC, 64 bit shift, and MWC
    :param seed: Seed to use
    :param m:
    :param a:
    :param c:
    :return:
    """

    # First linear congruential generator
    # While true, so the generator never stops making new numbers
    while True:
        # This is MLCC part
        generated_number = (a * seed + c) % m
        # Now bit shift
        generated_number = generated_number ^ (generated_number >> a1)
        generated_number = generated_number ^ (generated_number << a2)
        generated_number = generated_number ^ (generated_number >> a3)

        # Now MWC part
        mwc_out = a4 * (generated_number & (2 ** 32 - 1)) + (generated_number >> 32)

        seed = mwc_out
        # m = third_output
        mwc_out = mwc_out / m

        if mwc_out >= 1.:
            # Have to make it between 1 and 0, so mod 1. makes sure its between 0 and 1 now
            close_to_final = mwc_out % 1.
        else:
            close_to_final = mwc_out

        yield close_to_final


rand_gen = random_generator(seed)

first_thousand = []
first_thousand_x_1 = [0.0]
# 1b first one
for i in range(1000):
    first_thousand.append(next(rand_gen))
    print(first_thousand[i])
    if i > 0:
        first_thousand_x_1.append(first_thousand[i - 1])

# Now plot xi+1 vs xi

plt.scatter(first_thousand, first_thousand_x_1)
# plt.scatter([i for i in range(1000)], first_thousand_x_1)
plt.xlabel("$X_i$")
plt.ylabel("$X_{i+1}$")
plt.show()

first_million = []
for i in range(1000000):
    first_million.append(next(rand_gen))

plt.hist(first_million, bins=np.linspace(0.0, 1.0, 20))
plt.show()


def part_two(a, b, c, n=100, rand_generator=rand_gen):
    """
    Integrate density profile

    x = r/r_vir so th radius relative to the virial radius

    0 to x_max = 5 gives average total number of satellites


    n(x) = A<N_sat>(x/b)^(a-3) ep(-(x/b)^c)

    Solve lll_V n(r) dV = <N_sat>
    The 3D spherical integral from 0 to 5


    :param a:
    :param b:
    :param c:
    :param n:
    :return:
    """

    return NotImplementedError
