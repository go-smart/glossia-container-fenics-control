import dolfin as d

from problem import IREProblem

d.set_log_level(1000)

print "Initializing"

ire = IREProblem()

print "Solving"

ire.solve()

print "Plotting result"

ire.plot_result()
