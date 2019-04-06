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

seed = 5227

"""

Add x^2*4*pi for a 1D integral to make it the spherical one

Use midpoint rule, since there is a 0

A = 1/everything else


"""


def random_generator(seed, m=2 ** 64 - 1, a=2349543, c=913842, a1=21, a2=35, a3=4, a4=4294957665):
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
    piv = np.arange(0,n)
    for k in range(n-1):

        # piv
        max_row_index = np.argmax(abs(A[k:n,k])) + k
        piv[[k,max_row_index]] = piv[[max_row_index,k]]
        A[[k,max_row_index]] = A[[max_row_index,k]]

        # LU
        for i in range(k+1,n):
            A[i,k] = A[i,k]/A[k,k]
            for j in range(k+1,n):
                A[i,j] -= A[i,k]*A[k,j]

    return [A,piv]


def ufsub(L,b):
    """ Unit row oriented forward substitution """
    for i in range(L.shape[0]):
        for j in range(i):
            b[i] -= L[i,j]*b[j]
    return b


def bsub(U,y):
    """ Row oriented backward substitution """
    for i in range(U.shape[0]-1,-1,-1):
        for j in range(i+1, U.shape[1]):
            y[i] -= U[i,j]*y[j]
        y[i] = y[i]/U[i,i]
    return y

def bisect(arr, value):
    """
    Finds the index in the array closest to value
    :param arr:
    :param value:
    :return: Index of insertion point for the value in a sorted array
    """

    low = 0
    high = len(arr)
    while low < high:
        mid = (low+high)//2
        if value < arr[mid]:
            high = mid
        else:
            low = mid+1

    return low

def one_d_cube_spline(x, y):
    len_x = len(x)

    h = [x[i+1]-x[i] for i in range(len_x-1)]

    A = np.zeros((len_x,len_x))
    A[0,0] = 1.

    for i in range(len_x-1):
        if i != (len_x - 2):
            A[i+1,i+1] = 2*(h[i] + h[i+1]) # Going down the diagonal, do the 2*(differences)
        A[i+1, i] = h[i]
        A[i, i+1] = h[i] # The left and right sode of this one, which is the upper and lower edges

    A[0,1] = 0.0 # so natural cubic spline
    A[len_x - 1, len_x - 2] = 0.0
    A[len_x - 1, len_x - 1] = 1.0 # Cubic spline end, should be 0?

    B = np.zeros(len_x) # RHS of equation
    for i in range(len_x - 2):
        B[i + 1] = 3*(y[i+2] - y[i+1]) / h[i+1] - 3*(y[i+1] - y[i])/h[i]

    LU, piv = lu_factor(A)
    B = B[piv]
    ytmp = ufsub(LU, B)
    c = bsub(LU, ytmp)

    # Now can calculate B and D
    d = []
    b = []
    for i in range(len_x - 1):
        d.append((c[i + 1] - c[i]) / (3.0 * h[i]))
        tb = (y[i + 1] - y[i]) / h[i] - h[i] * \
             (c[i + 1] + 2.0 * c[i]) / 3.0
        b.append(tb)

    xs = np.arange(0, 5, 0.0001)
    print("X Min: {} Max: {}".format(x[0], x[-1]))
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

    plt.plot(xs, interpolated_points)
    plt.scatter(interp_data_points, measured_values, s=10)
    plt.xscale("log")
    plt.yscale("log")
    plt.show()
    plt.cla()

    return y, b, c, d

def estimate_with_spline(xs,y,b,c,d):
    interpolated_points = []
    import bisect
    for point in xs:
        point = np.log10(point)
        # Get closest point first
        if point < x[0]:
            interpolated_points.append(None)
            continue
        elif point > x[-1]:
            interpolated_points.append(None)
            continue
        i = bisect.bisect(x, point) - 1
        dx = point - x[i]
        interpolated_points.append(y[i] + b[i] * dx + c[i] * dx ** 2 + d[i] * dx ** 3)

x = [1e-4, 1e-2, 1e-1, 1, 5]
A = 1 / integration_alg(sat_equation_no_A, lower_bound=0, upper_bound=5, number_of_steps=10000)
y = [np.log10(sat_equation(r, A)) for r in x]

one_d_cube_spline(np.log10(x), y)
exit()


# TODO Add interpolation method here, cubic spline, etc.
def cubic_spline(xdata, ydata):
    """
    This method was chosen because in the slides it acheived fairly good fits hile being simpler than the Akima spline

    Natural because at i= and i = N-1, setting y'' = 0

    This is then the equation" y = Ayi + Byi+1 + Cy''i + Dy''i+1

    A = (xi+1 - x) / (xi + 1 - xi)
    B = 1 - A
    C = 1/6*(A^3-A)(xi+1-xi)^2
    D = 1/6*(B^3 - B)(xi+1-xi)^2

    first deriv = y2-y1/x2-x + (1-2t)*a(1-t)+bt/(x2-x1) + t(1-t)b-a/x2-x1

    second deriv 2*(b-2a+(a-b)*3t)/(x2-x1)^2)

    t(x) = x - x1/x2 - x1
    a =
    :return:
    """

    num_points = len(xdata) - 1

    # Assume sorted points, because these will be

    A = np.zeros(shape=(4*num_points,4*num_points))
    b = np.zeros(shape=(4*num_points,1)) # RHS

    for i in range(num_points):
        # 2n equations from
        A[i][4*i+0] = xdata[i]**3
        A[i][4*i+1] = xdata[i]**2
        A[i][4*i+2] = xdata[i]
        A[i][4*i+3] = 1
        b[i] = ydata[i]
        # This is x3+x2+x+1 = y

        A[num_points+i][4*i+0] = xdata[i+1]**3
        A[num_points+i][4*i+1] = xdata[i+1]**2
        A[num_points+i][4*i+2] = xdata[i+1]
        A[num_points+i][4*i+3] = 1
        b[num_points+i] = ydata[i+1]
        # This is xi+1**3+xi+1**2+xi+1+1 = yi+1

        if i == 0:
            continue

        # now for the other 2n-2 equations to make it cubic
        # 3*a*x**2 + 2*b*c + c = 3*ai+1*xi**2 + 2*bi+1*xi + ci+1
        # 6ai*xi + 2*bi = 6ai+1*xi + 2*bi+1
        # For 1 to n-1, so 2n-2 total
        A[2*num_points+(i-1)][4*(i-1)+0] = 3*xdata[i]**2
        A[2*num_points+(i-1)][4*(i-1)+1] = 2*xdata[i]
        A[2*num_points+(i-1)][4*(i-1)+2] = 1
        A[2*num_points+(i-1)][4*(i-1)+0+4] = -3*xdata[i]**2
        A[2*num_points+(i-1)][4*(i-1)+1+4] = -2*xdata[i]
        A[2*num_points+(i-1)][4*(i-1)+2+4] = -1
        b[2*num_points+(i-1)] = 0

        # Now for the second equation
        A[3*num_points+(i-1)][4*(i-1)+0] = 6*xdata[i]
        A[3*num_points+(i-1)][4*(i-1)+1] = 2
        A[3*num_points+(i-1)][4*(i-1)+0+4] = -6*xdata[i]
        A[3*num_points+(i-1)][4*(i-1)+1+4] = -2
        b[3*num_points+(i-1)] = 0

    # For natural spline, set endpoints
    A[3*num_points-1+0][0+0] += 6*xdata[0]
    A[3*num_points-1+0][0+1] += 2
    b[3*num_points-1+0] += 0

    A[3*num_points+num_points-1][4*(num_points-1)+0] += 6*xdata[num_points]
    A[3*num_points+num_points-1][4*(num_points-1)+1] += 2
    b[3*num_points+num_points-1] += 0

    from pprint import pprint
    pprint(A)
    print()
    print(b)

    import scipy.linalg
    x = scipy.linalg.solve(A, b)
    spline = []
    act_spline = []
    for i in range(num_points):
        spline.append({"u": xdata[i], "v": xdata[i+1],
                       "a": float(x[4*i+0]),
                       "b": float(x[4*i+1]),
                       "c": float(x[4*i+2]),
                       "d": float(x[4*i+3])})

    print(spline)

    for index, interval in enumerate(interp_data_points):
        if index == 0:
            continue
        for p in spline:
            if np.isclose(p['u'], interp_data_points[index-1]) and np.isclose(p['v'], interp_data_points[index]):
                xs = np.linspace(p['u'],p['v'])
                interpolated_points = [p['a']*i**3+p['b']*i**2+p['c']*i+p['d'] for i in xs]
                plt.plot(xs, interpolated_points)
                plt.scatter(interp_data_points, measured_values, s=20)
                plt.show()
                plt.cla()

    raise NotImplementedError

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
start_poisson = nums_in_bins[0]
for value in np.arange(start_poisson, 33):
    poisson_values.append(poisson(sum(nums_in_bins) / len(nums_in_bins), 0))
print("Poisson: {}".format(poisson_values))

#print(max(np.asarray(nums_in_bins)/sum(nums_in_bins)))
#print(min(np.asarray(nums_in_bins)/sum(nums_in_bins)))
#nums_in_bins = np.asarray(nums_in_bins)/sum(nums_in_bins)
bin_values, _, _ = plt.hist(nums_in_bins, bins=1000)
plt.cla()
bin_values = bin_values / sum(bin_values)
plt.hist(np.linspace(min(nums_in_bins), max(nums_in_bins), 1000), bins=1000, weights=bin_values)
plt.plot(np.arange(start_poisson, 33), poisson_values, 'r')
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
            indicies.append((i,j,k))

print(A_values.flatten().shape)

# TODO Save to file and load as it takes awhile to generate

# Now have the 3D cube of values, need to write an interpolator for them
# Can  take the 1D interpolator and apply it to here


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
