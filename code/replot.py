from __future__ import print_function

import numpy as N
import matplotlib.pyplot as P
import sys

import variables as v

if len(sys.argv) < 2:
    print("Must have a filename!")
    exit(1)

filename = sys.argv[1]

# Output cumulative coverage curves as CSV
output_arrays = N.loadtxt(filename)
output_arrays = zip(*output_arrays)

# Plot the coverage curves
for pair, a in zip(v.electrode_pairs, output_arrays[1:]):
    P.plot(output_arrays[0], a, label="%d - %d" % pair)

# Draw the plot
P.draw()
P.xlabel(r"Threshold level of $|E|$ ($\mathrm{J}$)")
P.xlim([0, 1000])
P.ylabel(r"Fraction of tumour beneath level")
P.ylim([0, 1])

# Show a legend for the plot
P.legend(loc=3)

# Display the plot
P.show(block=True)
