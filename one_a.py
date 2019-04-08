import numpy as np


def poisson(lam, k):
    """
    Possion Function (lambda^k)*(e^-lambda)/(k!)
    :param lam:
    :param k:
    :return:
    """
    # The exponential part, pretty much a direct translation from eqn. to code
    exponent_part = np.float64(np.exp(-1. * lam))
    first_part = np.float64(lam ** k)
    # Factorial part....
    fact = np.float64(1)
    for i in range(1, np.int64(k + 1)):
        fact *= i

    final_result = np.float64((first_part * exponent_part) / fact)
    return final_result

import sys
sys.stdout = open('1a.txt', 'w')
print("The Poisson value for $\lambda$ = {} and k = {} is: {} ".format(1, 0, poisson(1, 0)))
print("The Poisson value for $\lambda$ = {} and k = {} is: {} ".format(5, 10, poisson(5, 10)))
print("The Poisson value for $\lambda$ = {} and k = {} is: {} ".format(3, 21, poisson(2.6, 40)))
