import numpy as np


def two_c(A, a, b, c):
    import sys
    sys.stdout = open('2c.txt', 'w')

    def sat_equation(r, A, num_satallites=100):
        return A * num_satallites * (r / b) ** (a - 3) * np.exp(-(r / b) ** c)

    # Since we are integrating the sat_equation in 3D spherical integral, but only a dependence on r, integral
    # corresponds to sat_equation times dV, or the derivative of the volume of the sphere, so 4 * pi * r^2
    def three_d_integral(r, A, num_sats):
        volume_of_sphere = 4 * np.pi * r ** 2
        return volume_of_sphere * sat_equation(r, A, num_sats)

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
