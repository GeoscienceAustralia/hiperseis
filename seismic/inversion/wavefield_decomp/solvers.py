#!/usr/bin/env python
# coding: utf-8
"""
Objective function minimization solvers.
"""

class SolverGlobalMhMcmc:
    """
    Drop-in custom solver for scipy.optimize.minimize, based on Metrolpolis-Hastings Monte Carlo
    Markov Chain random walk with burn-in.

    Rather than returning one global solution, the solver returns the N best ranked local solutions.
    It also returns a probability distribution of each unknown based on Monte Carlo statistics.
    """
    pass
# end class
