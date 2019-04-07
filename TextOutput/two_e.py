import numpy as np
import matplotlib.pyplot as plt


def two_e(rand_gen, A, a, b, c):
    import sys
    sys.stdout = open('2e.txt', 'w')

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
    plt.ylabel("Log(Number of Satellites)")
    plt.xlabel("Log(X (R/R_vir))")
    plt.plot(np.arange(1e-4, 5, 0.001), n(np.arange(1e-4, 5, 0.001)), 'r', label='N(x)')
    plt.hist(np.logspace(-4, 0.69897000433, 20), weights=new_bin_values, bins=log_bins, label="Haloes")
    plt.legend(loc='best')
    plt.savefig("./plots/1000_haloes.png", dpi=300)
    plt.cla()

    return haloes, bin_values, log_bins
