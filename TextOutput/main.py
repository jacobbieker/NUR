import numpy as np
import matplotlib.pyplot as plt
from .one_a import poisson
from .one_b import random_generator, one_b
from .two_a import two_a
from .two_b import two_b
from .two_c import two_c
from .two_d import two_d
from .two_e import two_e
from .two_f import two_f
from .two_g import two_g

import sys
sys.stdout = open('seed.txt', 'w')


seed = 5227

print("The Seed for this project is: {}".format(seed))

rand_gen = random_generator(seed)

one_b(rand_gen)

A, a, b, c = two_a(rand_gen)

two_b(A, a, b, c)

two_c(A, a, b, c)

two_d(rand_gen, A, a, b, c)

haloes, bin_values, log_bins = two_e(rand_gen, A, a, b, c)

two_f(A, a, b, c)

two_g(haloes, bin_values, log_bins, poisson)

