import numpy as np

def two_a(rand_gen):
    import sys
    sys.stdout = open('2a.txt', 'w')
    # First part is generating a,b,c for the rest of part 2
    a = 1.4 * next(rand_gen) + 1.1
    b = 1.5 * next(rand_gen) + 0.5
    c = 2.5 * next(rand_gen) + 1.5

    print("The value for a is: {}".format(a))
    print("The value for b is: {}".format(b))
    print("The value for c is: {}".format(c))

    def sat_equation_no_A(r):
        """
        This is for finding A, only need the radius then
        :param r: Radius in terms of r = Radius/(Virial Radius)
        :return:
        """
        volume_of_sphere = 4 * np.pi * r ** 2
        return volume_of_sphere * ((r / b) ** (a - 3) * np.exp(-(r / b) ** c))


    def integration_alg(func, lower_bound, upper_bound, number_of_steps):
        """

        :param func: function to use in integration, when given a radius, will return the value at that point
        :param lower_bound: lower bound of integration
        :param upper_bound: upper bound of integration
        :param number_of_steps: number of steps to do
        :return:
        """

        # Current method if the midpoint rule, as the function is an improper integral from the lower bound being 0

        # Need to integrate from lower_bound to upper_bound in number_of_steps

        # lower bound is 0, if there is a radius of 0, then no satallites in it
        integration_value = 0
        step_size = (upper_bound - lower_bound) / number_of_steps  # The number of steps to take
        for i in range(number_of_steps):
            if i != 0:
                # Current step can be just i*step_size but only if integral always starts at 0
                # since it might not, need a starting point otherwise:
                current_step = lower_bound + i * step_size
                prev_step = lower_bound + (i - 1) * step_size

                # Current midpoint is the current step + the prev step divided by 2
                # F(mk) where mk = (tk + tk-1)/2
                current_midpoint = (current_step + prev_step) / 2
                integration_value += func(current_midpoint)

        # Last bit is the multiplication by the step size to get the full value
        integration_value *= step_size

        return integration_value


    # To get A, realize that <N_sat> is on both sides, so equation becomes
    # A * integral = 1, so A = 1 / integral, to do it then need different equation than in sat_equation
    A = 1 / integration_alg(sat_equation_no_A, lower_bound=0, upper_bound=5, number_of_steps=10000)
    print("A: ", A)

    return A, a,b,c
