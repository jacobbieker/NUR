import numpy as np
import matplotlib.pyplot as plt


def two_b(A,a,b,c):
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

    def sat_equation(r, A, num_satallites=100):
        return A * num_satallites * (r / b) ** (a - 3) * np.exp(-(r / b) ** c)


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
        plt.xscale("log")
        plt.yscale("log")
        plt.xlabel("Log(X (R/(Virial Radius))")
        plt.ylabel("Log(Number of Satellites)")
        plt.savefig("./plots/interpolation.png", dpi=300)
        plt.cla()

        return y, b, c, d

    x = [1e-4, 1e-2, 1e-1, 1, 5]
    A = 1 / integration_alg(sat_equation_no_A, lower_bound=0, upper_bound=5, number_of_steps=10000)
    y = [np.log10(sat_equation(r, A)) for r in x]

    y, b_interp, c_interp, d_interp = one_d_cube_spline(np.log10(x), y)
