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
    exponent_part = np.exp(-1.*lam)
    first_part = lam**k
     # Factorial part....
    fact = 1
    for i in range(1, k+1):
        fact *= i

    final_result = (first_part*exponent_part)/fact
    return final_result

print(poisson(1,0))
print(poisson(5,10))
print(poisson(3,20))
print(poisson(2.6,40))

seed = 1337


def random_generator(seed, m=2**64, a=2349543, c=913842, a1=21, a2=35, a3=4, a4=4294957665):
    """
    Generates psuedorandom numbers with a combination of (M)LCC, 64 bit shift, and MWC
    :param seed: Seed to use
    :param m:
    :param a:
    :param c:
    :return:
    """

    # First linear congruential generator, c = 0 -> MLCC
    while True:
        first_output = (a*seed + c) % m
        # Now bit shift
        second_output = first_output ^ (first_output >> a1)
        third_output = second_output ^ (second_output << a2)
        fourth_output = third_output ^ (third_output >> a3)

        # Now MWC part
        mwc_out = a4*(fourth_output & (2**32-1)) + (fourth_output >> 32)

        seed = mwc_out
        #m = third_output
        mwc_out = mwc_out / m

        if mwc_out >= 1.:
            # Have to make it between 1 and 0
            close_to_final = mwc_out % 1.
        else:
            close_to_final = mwc_out

        yield close_to_final


rand_gen = random_generator(seed)

first_thousand = []
first_thousand_x_1 = []
first_thousand_x_1.append(0.0)
# 1b first one
for i in range(1000):
    first_thousand.append(next(rand_gen))
    print(first_thousand[i])
    if i > 0:
        first_thousand_x_1.append(first_thousand[i-1])


# Now plot xi+1 vs xi

plt.scatter(first_thousand, first_thousand_x_1)
#plt.scatter([i for i in range(1000)], first_thousand_x_1)
plt.xlabel("$X_i$")
plt.ylabel("$X_{i+1}$")
plt.show()

first_million = []
for i in range(1000000):
    first_million.append(next(rand_gen))

plt.hist(first_million, bins=np.linspace(0.0,1.0,20))
plt.show()


def part_two(a,b,c,n=100, rand_generator=rand_gen):
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