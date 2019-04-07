import numpy as np
import matplotlib.pyplot as plt


def two_g(haloes, bin_values, log_bins, poisson):
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
            merge_sort(lower_half) # Recursively call until list is 1 element long
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
        haloes[index][0] = merge_sort(radiii) # Sort for easier time getting the elements

    selected_satellites = []

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
        selected_satellites.append(radiii[start_index:end_index])

    selected_satellites = np.asarray(selected_satellites)

    nums_in_bins = [len(sat) for sat in selected_satellites]
    nums_in_bins = merge_sort(nums_in_bins) # Sort so that median, and quartile are easy to calculate
    import sys
    sys.stdout = open('2g.txt', 'w')
    # From the slides
    print("Median Value: {}".format(nums_in_bins[int(len(nums_in_bins) / 2)]), flush=True)
    print("16th Percentile: {}".format(nums_in_bins[int(0.16 * len(nums_in_bins))]), flush=True)
    print("84th Percentile: {}".format(nums_in_bins[int(0.84 * len(nums_in_bins))]), flush=True)

    poisson_values = []
    start_poisson = nums_in_bins[-1]
    # Go from 0 to past the end to have the plot show a smooth function the whole plot
    for value in np.arange(0, start_poisson + 10):
        poisson_values.append(poisson(sum(nums_in_bins) / len(nums_in_bins), value))

    bins = np.arange(nums_in_bins[0], nums_in_bins[-1] + 1, 1)
    bin_values, _, _ = plt.hist(nums_in_bins, bins=bins)
    plt.cla()
    plt.xlim(nums_in_bins[0] - 1, nums_in_bins[-1] + 1)
    bin_values = bin_values / sum(bin_values)
    plt.hist(np.arange(nums_in_bins[0], nums_in_bins[-1], 1), bins=bins, weights=bin_values,
             label='Number of galaxies in bin')
    plt.plot(np.arange(0, start_poisson + 10), poisson_values, 'r', label='Poisson Distribution')
    plt.legend(loc='best')
    plt.ylabel("Probability")
    plt.xlabel("Number of galaxies in radial bin")
    plt.savefig("./plots/hist_poisson.png", dpi=300)
    plt.cla()
