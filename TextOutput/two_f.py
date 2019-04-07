import numpy as np


def two_f(A, a, b, c):
    import sys
    sys.stdout = open('2f.txt', 'w')

    def sat_equation(r, A, num_satallites=100):
        return A * num_satallites * (r / b) ** (a - 3) * np.exp(-(r / b) ** c)

    # Since we are integrating the sat_equation in 3D spherical integral, but only a dependence on r, integral
    # corresponds to sat_equation times dV, or the derivative of the volume of the sphere, so 4 * pi * r^2
    def three_d_integral(r, A, num_sats):
        volume_of_sphere = 4 * np.pi * r ** 2
        return volume_of_sphere * sat_equation(r, A, num_sats)

    def n(x):
        return three_d_integral(x, A, 100)

    def root_finder(bracket=[1e-8, 5], epsilon=0.001, max_iter=500, root_to_find=1 / 2):
        """
        Find the roots, easiest method to do is the bisection method, so deciding on that one

        First bracket the root, and see where the function changes sign -> that is the root

        Gaurunteed to work (according to the slides), unlike secant method or some of the other ones

        N(x) = y/2 => y = 2*N(x)

        Essentially finding two roots, need to bracket on one side with bracket of [0,Nmax], other with [Nmax,5]

        Uses root_to_find by subtracting that value from the rest of the stuff, so that the root wanted is at 0

        :return:
        """

        for i in range(max_iter):
            mid_point = (bracket[0] + bracket[1]) / 2.
            # Now check to see which half the root is in
            # Is the midpoint the root?
            mid_value = n(mid_point) - root_to_find
            f_a = n(bracket[0]) - root_to_find
            if f_a * mid_value < 0:  # They have opposite signs
                bracket[1] = mid_point
            elif -1. * epsilon < mid_value < epsilon:  # This is the root within epsilon
                return mid_point
            else:
                # Other two have opposite signs
                bracket[0] = mid_point
            # Now check if it has converged, that is if the space between the brackets is within epsilon
        # If it gets to here, then no root found after the max iterations
        print("Root not found in {} iterations".format(max_iter))
        return (bracket[0] + bracket[1]) / 2.

    Nmax = n(np.arange(1e-4, 5, 0.001))
    max_index = list(Nmax).index(max(Nmax))
    max_val = np.arange(1e-4, 5, 0.001)[max_index]

    # There are two roots in this problem, as half the max happens on either side of the max
    lower_root = root_finder(bracket=[1e-8, max_val], root_to_find=n(max_val) / 2.)
    upper_root = root_finder(bracket=[max_val, 5], root_to_find=n(max_val) / 2.)
    print("Roots = {}, {}".format(lower_root, upper_root))
