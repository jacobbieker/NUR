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
# print(poisson(101,200))

seed = 1337

"""

Add x^2*4*pi for a 1D integral to make it the spherical one

Use midpoint rule, since there is a 0

A = 1/everything else


"""


def random_generator(seed, m=2 ** 64, a=2349543, c=913842, a1=21, a2=35, a3=4, a4=4294957665):
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

"""

End of Part 1

"""


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

    def interpolation_aspect():
        """
        This performs the interpolation part of 2 b)
        :return:
        """
        r_val = [1E-4, 1E-2, 1E-1, 1, 5]
        integral_vals = [sat_equation_no_A(r) for r in r_val]

    def derivative_dn_dx(b):
        """
        Returns the numerical derivative of the equation at x = b
        Returns the analytical result as well
        :param b:
        :return:
        """

        raise NotImplementedError

    def generate_satallite_positions(number_of_satallites):
        """
        Generates random positions that follow the profile of equation 2

        To sample the distribution, rejection sampling is used. The reason for this is the ease of implementing it
        for this problem, since only have to check if the random sample is less than the p(x) given in the handin

        In rejetion sampling, you accept the x value if the y value for that x is less than or equal to p(x)

        p(x) in this case is n(x)*4*pi*x^2 dx / N_sat = (4*pi*A*b^3*(x/b)^a*exp(-(x/b)^c)*(a-c*(x/b)^c-1)/(x^2))
        The N_sat cancels with the one in the n(x)

        :param number_of_satallites: Number of satallites to generate
        :return:
        """

        def p_x(x):
            """
            Since only using 100 satallites for each, this is the way to go
            :param x:
            :return:
            """
            return three_d_integral(x, A, num_sats=100)

        # Now generate random x values, plug into sat equation to get a value, and see if its less than p_x
        sat_radius = []
        sat_input = []
        sat_theta = []
        sat_omega = []

        # Need min and max values for the rejection y sample
        min_y_value = min(
            [sat_equation(i, A, num_satallites=number_of_satallites) / number_of_satallites for i in r_val])

        while len(sat_radius) < number_of_satallites:
            input_x = next(rand_gen) * 5
            random_radius = sat_equation(input_x, A=A, num_satallites=number_of_satallites) / number_of_satallites
            if random_radius <= p_x(random_radius):
                sat_radius.append(random_radius)
                sat_input.append(input_x)


"""


Part 2. This part deals with a 3D spherical integral and other things


"""

# First part is generating a,b,c for the rest of part 2
a = 1.4 * next(rand_gen) + 1.1
b = 1.5 * next(rand_gen) + 0.5
c = 2.5 * next(rand_gen) + 1.5


def sat_equation(r, A, num_satallites=100):
    return A * num_satallites * (r / b) ** (a - 3) * np.exp(-(r / b) ** c)

    # Since we are integrating the sat_equation in 3D spherical integral, but only a dependence on r, integral
    # corresponds to sat_equation times dV, or the derivative of the volume of the sphere, so 4 * pi * r^2


def three_d_integral(r, A, num_sats):
    volume_of_sphere = 4 * np.pi * r ** 2
    return volume_of_sphere * sat_equation(r, A, num_sats)


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
    step_size = (upper_bound - lower_bound) / number_of_steps
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

interp_data_points = [1e-4, 1e-2, 1e-1, 1, 5]

measured_values = [sat_equation(r, A) for r in interp_data_points]


# Now need to interpolate between the points

# TODO Add interpolation method here, cubic spline, etc.

# Now onto Part c, numerical differentiation

def derivative(func, b, step_size=0.001):
    """
    This uses the central differences method to calculate the derivative of a function
    :param b:
    :return:
    """

    deriv = func(b + step_size) - func(b - step_size)
    deriv = deriv / (2 * step_size)
    return deriv


def analytic_derivative(b):
    """
    Computed on Wolfram Alpha, this is the actual derivative with respeict to dx

    Since at x =b ,
    :return:
    """
    x = b
    return (4 * np.pi * A * b ^ 3 * (x / b) ^ a * np.exp(-(x / b) ^ c) * (a - c * (x / b) ^ c - 1) / (x ^ 2))


# TODO Compare the Two for X = b

# Part D Sampling

def random_sample(func, xmin, xmax, ymin, ymax, num_samples):
    """
    This random sampling uses the rejection method, primarily because its the easiest to implement

    For this project, since x can be between 0 and 5, y is also between 0 and 5

    Generate random numbers in both, if y < p_x(x) then the data point is accepted

    :return:
    """

    inputs = []
    outputs = []

    while len(outputs) < num_samples:
        x = next(rand_gen) * (xmax - xmin) + xmin
        y = next(rand_gen) * (ymax - ymin) + ymin
        if y <= func(x, A, 1):
            inputs.append(x)
            outputs.append(y)

    return inputs, outputs


rand_sample_x, rand_sample_y = random_sample(three_d_integral, 0, 5, 0, 5, 10000)
print("Max: ", np.max(rand_sample_x))
plt.scatter(rand_sample_x, rand_sample_y, s=1)
plt.plot(np.arange(0, 5, 0.001), [three_d_integral(i, A, 1) for i in np.arange(0, 5, 0.001)], 'r')
plt.title("Random Sampling")
plt.show()


# Now need to create random phi and theta for the given r values

def create_halo(number_of_satallites):
    rand_sample_x = random_sample(three_d_integral, 0, 5, 0, 5, number_of_satallites)[0]

    # now randomize the phi and thetas
    phi_sample = []
    theta_sample = []

    for _ in rand_sample_x:
        phi_sample.append((2 * np.pi * next(rand_gen)))
        theta_sample.append(np.pi * next(rand_gen))

    return rand_sample_x, phi_sample, theta_sample


# Part e, the log-log histogram

# TODO Make log log histogram

# Given N(x)

def n(x):
    return three_d_integral(x, A, 100)

plt.title("Log-Log of N(x)")
plt.xscale("log")
plt.yscale("log")
plt.plot(np.arange(1e-4,5,0.001), n(np.arange(1e-4,5,0.001)), 'r')
plt.show()

# Need to Generate 1000 Halos Now

# TODO Generate 100 halos

"""

Part f Root Finding

"""


def root_finder():
    """
    Find the roots, easiest method to do is the bisection method, so deciding on that one

    :return:
    """



    raise NotImplementedError



def cubic_spline():
    raise NotImplementedError


def central_difference():
    raise NotImplementedError


def analyitic_N(x):
    raise NotImplementedError




def p_x(x):
    raise NotImplementedError


def interpolation():
    raise NotImplementedError
