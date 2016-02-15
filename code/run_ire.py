from __future__ import print_function

import fenics as d

print("Importing...")

from problem import IREProblem

d.set_log_level(1000)

print("Initializing")

ire = IREProblem()

print("Loading")

ire.load()

print("Solving")

ire.solve()

print("Plotting result")

#ire.plot_bitmap_result()
ire.plot_result()
