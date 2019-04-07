import numpy as np


def two_h():
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

    # TODO Save to file and load as it takes awhile to generate

    # Now have the 3D cube of values, need to write an interpolator for them
    # Can  take the 1D interpolator and apply it to here

    # Have to do it N times, where N is the number of AxB for C values, so huge
    # Unless take the closest 8 points
    subcube = A_values[int(len(a_range) / 2 - 2):int(len(a_range) / 2 + 2),
              int(len(b_range) / 2 - 2):int(len(b_range) / 2 + 2), int(len(b_range) / 2 - 2):int(len(b_range) / 2 + 2)]


    # So now try interpolating points based off the 64 points around it, allows for
    # 3D interpolation = 3 1D interpolations
    # Or 3D linear interpolations


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
