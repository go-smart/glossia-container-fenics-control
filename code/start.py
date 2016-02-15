from __future__ import print_function

from gosmart import setup as setup_gssa
import os
setup_gssa(os.path.join(os.path.dirname(__file__), 'parameters.yml'))

from go_smart_ire.problem import IREProblem
import fenics as d


def main():

    d.set_log_level(1000)

    print("Initializing")

    ire = IREProblem()

    print("Loading")

    ire.load()

    print("Solving")

    ire.solve()

    #print("Plotting result")

    #ire.plot_result()

    print("Saving output lesion")

    ire.save_lesion()

if __name__ == "__main__":
    main()
