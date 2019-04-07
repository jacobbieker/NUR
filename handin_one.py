import numpy as np
import matplotlib.pyplot as plt

seed = 5227

print("The Seed for this project is: {}".format(seed))


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


print("The Poisson value for $\lambda$ = {} and k = {} is: {} ".format(1, 0, poisson(1, 0)))
print("The Poisson value for $\lambda$ = {} and k = {} is: {} ".format(5, 10, poisson(5, 10)))
print("The Poisson value for $\lambda$ = {} and k = {} is: {} ".format(3, 21, poisson(2.6, 40)))


# print(poisson(101,200))


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


rand_gen = random_generator(seed)

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
plt.show()

first_million = []
for i in range(1000000):
    first_million.append(next(rand_gen))

plt.hist(first_million, bins=np.linspace(0.0, 1.0, 20))
plt.show()

"""

End of Part 1

"""

"""


Part 2. This part deals with a 3D spherical integral and other things


"""

# First part is generating a,b,c for the rest of part 2
a = 1.4 * next(rand_gen) + 1.1
b = 1.5 * next(rand_gen) + 0.5
c = 2.5 * next(rand_gen) + 1.5

print("The value for a is: {}".format(a))
print("The value for b is: {}".format(b))
print("The value for c is: {}".format(c))


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

interp_data_points = [1e-4, 1e-2, 1e-1, 1, 5]

measured_values = [np.log10(sat_equation(r, A)) for r in interp_data_points]


# Now need to interpolate between the points
def lu_factor(A):
    """
        LU factorization with partial pivoting

        Overwrite A with:
            U (upper triangular) and (unit Lower triangular) L
        Return [LU,piv]
            Where piv is 1d numpy array with row swap indices
    """
    n = A.shape[0]
    piv = np.arange(0, n)
    for k in range(n - 1):

        # piv
        max_row_index = np.argmax(abs(A[k:n, k])) + k
        piv[[k, max_row_index]] = piv[[max_row_index, k]]
        A[[k, max_row_index]] = A[[max_row_index, k]]

        # LU
        for i in range(k + 1, n):
            A[i, k] = A[i, k] / A[k, k]
            for j in range(k + 1, n):
                A[i, j] -= A[i, k] * A[k, j]

    return [A, piv]


def ufsub(L, b):
    """ Unit row oriented forward substitution """
    for i in range(L.shape[0]):
        for j in range(i):
            b[i] -= L[i, j] * b[j]
    return b


def bsub(U, y):
    """ Row oriented backward substitution """
    for i in range(U.shape[0] - 1, -1, -1):
        for j in range(i + 1, U.shape[1]):
            y[i] -= U[i, j] * y[j]
        y[i] = y[i] / U[i, i]
    return y


def bisect(arr, value):
    """
    Finds the index in the array closest to value
    :param arr: The array of values
    :param value: The value to insert/the interpolation value
    :return: Index of insertion point for the value in a sorted array
    """

    low = 0
    high = len(arr)
    while low < high:
        mid = int((low + high) / 2)  # Get the midpoint to test if the value is above or below it
        if value < arr[mid]:
            high = mid
        else:
            low = mid + 1
    return low


def one_d_cube_spline(x, y):
    """
    This method was chosen because in the slides it achieved fairly good fits while being simpler than the Akima spline

    Natural because at i= and i = N-1, setting y'' = 0

    This is then the equation" y = Ayi + Byi+1 + Cy''i + Dy''i+1

    A = (xi+1 - x) / (xi + 1 - xi)
    B = 1 - A
    C = 1/6*(A^3-A)(xi+1-xi)^2
    D = 1/6*(B^3 - B)(xi+1-xi)^2

    first deriv = y2-y1/x2-x + (1-2t)*a(1-t)+bt/(x2-x1) + t(1-t)b-a/x2-x1

    second deriv 2*(b-2a+(a-b)*3t)/(x2-x1)^2)

    t(x) = x - x1/x2 - x1
    :return:
    """
    len_x = len(x)

    h = [x[i + 1] - x[i] for i in range(len_x - 1)]

    A = np.zeros((len_x, len_x))
    A[0, 0] = 1.

    for i in range(len_x - 1):
        if i != (len_x - 2):
            A[i + 1, i + 1] = 2 * (h[i] + h[i + 1])  # Going down the diagonal, do the 2*(differences)
        A[i + 1, i] = h[i]
        A[i, i + 1] = h[i]  # The left and right sode of this one, which is the upper and lower edges

    A[0, 1] = 0.0  # so natural cubic spline
    A[len_x - 1, len_x - 2] = 0.0
    A[len_x - 1, len_x - 1] = 1.0  # Cubic spline end, should be 0?

    B = np.zeros(len_x)  # RHS of equation
    for i in range(len_x - 2):
        # This is the
        B[i + 1] = 3 * (y[i + 2] - y[i + 1]) / h[i + 1] - 3 * (y[i + 1] - y[i]) / h[i]

    LU, piv = lu_factor(A) # Get the LU decomposition and pivot indicies
    B = B[piv] # reorder B to match the pivots
    ytmp = ufsub(LU, B) # Do forward substitution to get y for y = ax
    c = bsub(LU, ytmp) # Do backward substitution to get x from x = LU*y

    # Now can calculate B and D
    d = []
    b = []
    for i in range(len_x - 1):
        d.append((c[i + 1] - c[i]) / (3.0 * h[i]))
        tb = (y[i + 1] - y[i]) / h[i] - h[i] * \
             (c[i + 1] + 2.0 * c[i]) / 3.0
        b.append(tb)

    xs = np.arange(0, 5, 0.0001)
    interpolated_points = []
    for point in xs:
        point = np.log10(point)
        # Get closest point first
        if point < x[0]:
            interpolated_points.append(None)
            continue
        elif point > x[-1]:
            interpolated_points.append(None)
            continue
        i = bisect(x, point) - 1 # Find the closest point through determining where the input falls in the array
        dx = point - x[i] # Difference between the measured point and the interpolation point
        interpolated_points.append(y[i] + b[i] * dx + c[i] * dx ** 2 + d[i] * dx ** 3)
        # Uses the coefficients to create y + b*x + c*x^2 + d*x^3 = x

    plt.plot(xs, interpolated_points, c='r', label='Interpolated Values')
    plt.scatter(interp_data_points, measured_values, s=10, label='Measured Values')
    plt.legend(loc='best')
    plt.xlabel("X (R/(Virial Radius)")
    plt.ylabel("Number of Satellites")
    plt.title("")
    plt.show()
    plt.cla()

    return y, b, c, d


def estimate_with_spline(xs, y, b, c, d):
    interpolated_points = []
    for point in xs:
        point = np.log10(point)
        # Get closest point first
        if point < x[0]:
            interpolated_points.append(None)
            continue
        elif point > x[-1]:
            interpolated_points.append(None)
            continue
        i = bisect(x, point) - 1
        dx = point - x[i]
        interpolated_points.append(y[i] + b[i] * dx + c[i] * dx ** 2 + d[i] * dx ** 3)

    return xs, interpolated_points


x = [1e-4, 1e-2, 1e-1, 1, 5]
A = 1 / integration_alg(sat_equation_no_A, lower_bound=0, upper_bound=5, number_of_steps=10000)
y = [np.log10(sat_equation(r, A)) for r in x]

y, b_interp, c_interp, d_interp = one_d_cube_spline(np.log10(x), y)
xs, interpolated_points = estimate_with_spline(xs=np.arange(0, 5, 0.0001), y=y, b=b_interp, c=c_interp, d=d_interp)

plt.plot(xs, interpolated_points)
plt.scatter(interp_data_points, measured_values, s=10)
# plt.xscale("log")
# plt.yscale("log")
plt.show()
plt.cla()


# Now onto Part c, numerical differentiation

def n(x):
    return three_d_integral(x, A, 1)


def analytic_derivative(b):
    """
    Computed on Wolfram Alpha, this is the actual derivative with respeict to dx

    Since at x =b ,
    :return:
    """
    x = b
    return (4 * np.pi * A * b ** 3 * (x / b) ** a * np.exp(-(x / b) ** c) * (a - c * (x / b) ** c - 1) / (x ** 2))


def derivative(func, b, step_size=0.1, iterations=5):
    """
    This uses the central differences method to calculate the derivative of a function

    The step size was chosen to minimize the error between the numerical and analytical results, smaller step size resulted
    in a larger error, as well as a larger step size

    Ridder method: keep decreasing step_size until error grows

    A(1,m) = f(x+h/2^m-1) - f(x-h/2^m-1)/(2h/s^m-1)
    A(n,m) = 4^n-1*A(n-1,m+1) - A(n-1,m)/(4^n-1)-1)

    :param b:
    :return:
    """

    def A_deriv(n, m):
        if n == 1:
            result = (func(b + step_size / 2 ** (m - 1)) - func(b - step_size / (2 ** (m - 1)))) / (
                    2 * step_size / (2 ** (m - 1)))
        else:
            result = (4 ** (n - 1) * A_deriv(n - 1, m + 1) - A_deriv(n - 1, m)) / (4 ** (n - 1) - 1)
        return result

    best_approx = np.inf
    m = 1
    for i in range(iterations):
        deriv = A_deriv(i + 1, m)
        m += 1
        if analytic_derivative(b) - deriv < best_approx:
            best_approx = deriv

    return deriv


print("Analytic: {}\n Numerical: {}\n Difference: {}\n".format(np.round(analytic_derivative(b), 12),
                                                               np.round(derivative(n, b), 12),
                                                               np.round(analytic_derivative(b), 12) - np.round(
                                                                   derivative(n, b), 12)))


# Part D Sampling

def random_sample(func, xmin, xmax, ymin, ymax, num_samples):
    """

            Generates random positions that follow the profile of equation 2

        To sample the distribution, rejection sampling is used. The reason for this is the ease of implementing it
        for this problem, since only have to check if the random sample is less than the p(x) given in the handin

        In rejetion sampling, you accept the x value if the y value for that x is less than or equal to p(x)

        p(x) in this case is n(x)*4*pi*x^2 dx / N_sat = (4*pi*A*b^3*(x/b)^a*exp(-(x/b)^c)*(a-c*(x/b)^c-1)/(x^2))
        The N_sat cancels with the one in the n(x)

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


rand_sample_x, rand_sample_y = random_sample(three_d_integral, 1e-8, 5, 1e-8, 5, 10000)
print("Max: ", np.max(rand_sample_x))
plt.scatter(rand_sample_x, rand_sample_y, s=1)
plt.plot(np.arange(0, 5, 0.001), [three_d_integral(i, A, 1) for i in np.arange(0, 5, 0.001)], 'r')
plt.title("Random Sampling")
# plt.xscale('log')
# plt.yscale('log')
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


# Now outputting them
# for r, p, t in create_halo(100):
#    print("R: {} $\phi$: {} $\\theta$: {}".format(r, p, t))


# Part e, the log-log histogram


# Given N(x)

def n(x):
    return three_d_integral(x, A, 100)


# Need to Generate 1000 Halos Now

log_bins = np.logspace(-4, 0.69897000433, 21)  # Need 21 for the end of the bins, goes from 1e-4 to 5 in logspace


def create_haloes(number_of_haloes):
    """
    Creates a set number of haloes with 100 satallites each
    :param number_of_haloes:
    :return:
    """

    haloes = []
    radii = []
    for i in range(number_of_haloes):
        r, p, t = create_halo(100)
        haloes.append([r, p, t])
        radii.append(r)

    radii = np.asarray(radii)
    radii = np.concatenate(radii)
    return haloes, radii


def calc_avg_satallites_per_bin(bin_values, bins, num_haloes):
    """
    Divide bin values by width of bins, then by number of haloes used to create it
    Gives average number of satallies per bin
    :param bin_values:
    :param bins:
    :return:
    """

    new_averages = []

    for index, element in enumerate(bin_values):
        avg_per_halo = element / num_haloes  # Divide by number of haloes to get average per halo
        avg_per_bin_width = avg_per_halo / (
                bins[index + 1] - bins[index])  # Divide by bin width to normalize for bin width
        new_averages.append(avg_per_bin_width)

    return np.asarray(new_averages)


haloes, radii = create_haloes(1000)
bin_values, bins, _ = plt.hist(radii, bins=log_bins)
plt.cla()
new_bin_values = calc_avg_satallites_per_bin(bin_values, bins, 1000)
plt.title("Log-Log of N(x)")
plt.xscale("log")
plt.yscale("log")
plt.plot(np.arange(1e-4, 5, 0.001), n(np.arange(1e-4, 5, 0.001)), 'r')
plt.hist(np.logspace(-4, 0.69897000433, 20), weights=new_bin_values, bins=log_bins)
plt.show()

"""

Part f Root Finding

"""


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
print("Root = {}, {}".format(lower_root, upper_root))
print("N(x) Value at Roots: {} {} \n Max Value: {}".format(n(lower_root), n(upper_root), n(max_val)))

"""

Part g) Sorting, histogram, and Poisson checking

"""

# Get the radial bin with the largest number of galaxies
index_of_radial_bin_max = list(bin_values).index(max(bin_values))
inner_radius = log_bins[index_of_radial_bin_max]
outer_radius = log_bins[index_of_radial_bin_max + 1]


# Now sort the 1000 haloes on their ownm using Mergesort
def merge_sort(arr):
    """

    Mergesort implementation

    Takes an array, which is then split recursively until there is only single element arrays

    Then builds them back up, sorting them as it goes, resulting in O(nlogn) time
    """

    if len(arr) > 1:
        # Split array
        mid_point = int(len(arr) / 2.)
        lower_half = arr[:mid_point]
        upper_half = arr[mid_point:]
        merge_sort(lower_half)
        merge_sort(upper_half)

        # This only occurs after arr is split into all 1 element arrays
        # Starts with 2 subarrays of 1, building up to the whole list

        i = 0
        j = 0
        k = 0
        while i < len(lower_half) and j < len(upper_half):
            if lower_half[i] < upper_half[j]:
                arr[k] = lower_half[i]
                i += 1
            else:
                arr[k] = upper_half[j]
                j += 1
            k += 1

        # Now if N is not even, then either the lower or upper half will have extra elements
        # so need to add those, already in order, elements
        for l in range(i, len(lower_half)):
            arr[k] = lower_half[l]
            k += 1
        for l in range(j, len(upper_half)):
            arr[k] = upper_half[l]
            k += 1

        # Now arr is sorted, return it
        return arr


# Then only select ones that fall within that range between the two
for index, elements in enumerate(haloes):
    radiii, _, _ = elements
    haloes[index][0] = merge_sort(radiii)

selected_satallites = []

for radiii, _, _ in haloes:
    # Now they are sorted, only select the ones inbetween the radial values
    start_index = -1
    end_index = -1
    for index, element in enumerate(radiii):
        if element > inner_radius and start_index < 0:
            start_index = index
            break
    # Reversing list and going backwards means we can ignore where most of the galaxies are, and only select the one
    # that creates the border of the bin
    for index, element in reversed(list(enumerate(radiii))):
        if element < outer_radius and end_index < 0:
            end_index = index + 1  # Need to add one since slice does not include the last element specified
            break
    selected_satallites.append(radiii[start_index:end_index])

selected_satallites = np.asarray(selected_satallites)

nums_in_bins = [len(sat) for sat in selected_satallites]
nums_in_bins = merge_sort(nums_in_bins)

print("Median Value: {}".format(nums_in_bins[int(len(nums_in_bins) / 2)]))
print("16th Percentile: {}".format(nums_in_bins[int(0.16 * len(nums_in_bins))]))
print("84th Percentile: {}".format(nums_in_bins[int(0.84 * len(nums_in_bins))]))

poisson_values = []
start_poisson = nums_in_bins[-1]
for value in np.arange(0, start_poisson + 10):
    poisson_values.append(poisson(sum(nums_in_bins) / len(nums_in_bins), value))
print("Poisson: {}".format(poisson_values))

# print(max(np.asarray(nums_in_bins)/sum(nums_in_bins)))
# print(min(np.asarray(nums_in_bins)/sum(nums_in_bins)))
# nums_in_bins = np.asarray(nums_in_bins)/sum(nums_in_bins)
bins = np.arange(min(nums_in_bins), max(nums_in_bins) + 1, 1)
bin_values, _, _ = plt.hist(nums_in_bins, bins=bins)
plt.cla()
plt.xlim(min(nums_in_bins) - 1, max(nums_in_bins) + 1)
bin_values = bin_values / sum(bin_values)
plt.hist(np.arange(min(nums_in_bins), max(nums_in_bins), 1), bins=bins, weights=bin_values)
plt.plot(np.arange(0, start_poisson + 10), poisson_values, 'r')
plt.show()

"""

Part 2h interpolate in 3D

"""


def sat_equation_A_cube(r, a, b, c):
    """
    This is for finding A, only need the radius then
    :param r: Radius in terms of r = Radius/(Virial Radius)
    :return:
    """
    volume_of_sphere = 4 * np.pi * r ** 2
    return volume_of_sphere * ((r / b) ** (a - 3) * np.exp(-(r / b) ** c))


def integration_alg(func, lower_bound, upper_bound, number_of_steps, a, b, c):
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
            integration_value += func(current_midpoint, a, b, c)

    # Last bit is the multiplication by the step size to get the full value
    integration_value *= step_size

    return integration_value

    # To get A, realize that <N_sat> is on both sides, so equation becomes
    # A * integral = 1, so A = 1 / integral, to do it then need different equation than in sat_equation


# First need to create the 3D cube of values
b_range = np.arange(0.5, 2.1, 0.1)
a_range = np.arange(1.1, 2.6, 0.1)
c_range = np.arange(1.5, 4.1, 0.1)

A_values = np.zeros((len(a_range), len(b_range), len(c_range)))
indicies = []
# now fill in the 6240 values in the 3D cube of values
for i, a in enumerate(a_range):
    for j, b in enumerate(b_range):
        for k, c in enumerate(c_range):
            A_values[i, j, k] = 1 / integration_alg(sat_equation_A_cube,
                                                    lower_bound=0,
                                                    upper_bound=5,
                                                    number_of_steps=10000,
                                                    a=a, b=b, c=c)
            indicies.append((i, j, k))

print(A_values.flatten().shape)

# TODO Save to file and load as it takes awhile to generate

# Now have the 3D cube of values, need to write an interpolator for them
# Can  take the 1D interpolator and apply it to here

# Have to do it N times, where N is the number of AxB for C values, so huge
# Unless take the closest 8 points
subcube = A_values[int(len(a_range) / 2 - 2):int(len(a_range) / 2 + 2),
          int(len(b_range) / 2 - 2):int(len(b_range) / 2 + 2), int(len(b_range) / 2 - 2):int(len(b_range) / 2 + 2)]

print(subcube.shape)
print(subcube.flatten().shape)
print(subcube)


# So now try interpolating points based off the 64 points around it, allows for
# 3D interpolation = 3 1D interpolations
# Or 3D linear interpolations

def interpolator_3d(a, b, c, cube, size_subcube=2):
    """

    Interpolator that interpolates in 3D when given a data cube and an x,y,z point in a,b,c range

    :param cube:
    :return:
    """
    # To handle various resolutions of cubes
    a_range = np.linspace(1.1, 2.6, cube.shape[0])
    b_range = np.linspace(0.5, 2.1, cube.shape[1])
    c_range = np.linspace(1.5, 4.1, cube.shape[2])

    # Find the value, using bisect in each dimension
    a_loc = bisect(a_range, a) - 1
    b_loc = bisect(b_range, b) - 1
    c_loc = bisect(c_range, c) - 1

    # now gather the subcube around it
    # First need to add another layer outside the cube for edge cases. Because we assume its a straight line outside the
    # spline, we can just copy all the values out one more
    # Makes sure that there is no index out of bounds for this
    cube = np.pad(cube, pad_width=size_subcube, mode="edge")
    subcube = cube[a_loc - size_subcube:a_loc + size_subcube, b_loc - size_subcube:b_loc + size_subcube,
              c_loc - size_subcube:c_loc + size_subcube]

    # Now onto the actual spline interpolation in 3D
    # TODO Add 3D spline interpolation based on the subcube
    # Because there is a spline for each single-width vector in the cube, the number of splines goes up as
    # N^2, because, for example, adding a single more a column in the interpolation means that b*c more splines must be
    # generated
    # (1) Per-formMspline interpolations to get a vector of valuesy.x1i;x2/,iD0;:::;M1.
    # (2)  Construct  a  one-dimensional  spline  through those  values.
    # (3)  Finally,  spline-interpolate to the desired valuey.x1;x2/

    # So plan is to do M splines through c_range first, to get a 2D array of y(ai,bi,c)
    # Then N splines through b space to get 1D array of y(ai,b,c)
    # Finally, 1D spline through ai to get the final value of y(a,b,c) = A


"""

Part 3

"""

# Read in file, one at a time for memory limits

with open("satgals_m12.txt", "r") as mass_haloes:
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

"""

b Fit a function. because otherwise interpolating 3 values based off a single value, whle fitting a function to each
should be easier to do.

"""
