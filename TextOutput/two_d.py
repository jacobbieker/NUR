import numpy as np
import matplotlib.pyplot as plt


def two_d(rand_gen, A, a, b, c):
    import sys
    sys.stdout = open('2d.txt', 'w')

    def sat_equation(r, A, num_satallites=100):
        return A * num_satallites * (r / b) ** (a - 3) * np.exp(-(r / b) ** c)

    # Since we are integrating the sat_equation in 3D spherical integral, but only a dependence on r, integral
    # corresponds to sat_equation times dV, or the derivative of the volume of the sphere, so 4 * pi * r^2
    def three_d_integral(r, A, num_sats):
        volume_of_sphere = 4 * np.pi * r ** 2
        return volume_of_sphere * sat_equation(r, A, num_sats)

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
    plt.scatter(rand_sample_x, rand_sample_y, s=1, label='Sampled Points')
    plt.plot(np.arange(0, 5, 0.001), [three_d_integral(i, A, 1) for i in np.arange(0, 5, 0.001)], 'r', label='p(x)')
    plt.legend(loc='best')
    plt.title("Random Sampling")
    plt.xlabel("X (R/R_vir)")
    plt.ylabel("p(x)")
    plt.savefig("./plots/random_sample.png", dpi=300)
    plt.cla()

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
    print("(R, $\phi$, $\\theta$)")
    for r, p, t in create_halo(100):
        print("{}, {}, {}".format(r, p, t))
