# Script with functions to find values for a and b in the new TT method
from scipy.optimize import root, least_squares
from random import uniform
#from sklearn.cluster import Birch
from time import time
import numpy as np


def equations(p):
    a, b = p
    return (a*b*(a+1)/((a+b)*(a+b+1)*(a+b+2))-r,
            a*b*(a+1)*(b+1)/((a+b)*(a+b+1)*(a+b+2)*(a+b+3))-s)


def r_s_equation(a,b):
    return (a*b*(a+1)/((a+b)*(a+b+1)*(a+b+2)),
            a*b*(a+1)*(b+1)/((a+b)*(a+b+1)*(a+b+2)*(a+b+3)))


def find_roots(r_new, s_new):
    global r, s
    r = r_new
    s = s_new
    roots = []
    eq_values = []
    tol = 1e-9
    search_space = [[0, 0], [1, 10000]]
    i=0
    while i<500 or len(roots)==0:
        i = i + 1
        start_guess = [uniform(search_space[0][0], search_space[1][0]),
                       uniform(search_space[0][1], search_space[1][1])]
        a, b = root(equations, start_guess, method='hybr', tol=1e-9).x
        eq_1, eq_2 = equations((a, b))
        eq_tot = abs(eq_1)+abs(eq_2)
        if (eq_tot) < tol and [a, b] not in roots:
            roots.append((a, b))
            eq_values.append(eq_tot)
        if i > 1000:
            #raise Exception("Root values could not be found within reasonable time.")
            roots.append((np.nan, np.nan))
            eq_values.append(np.nan)

    #brc = Birch(n_clusters=None)
    #brc.fit(roots)
    #return brc.subcluster_centers_

    index = eq_values.index(min(eq_values))
    return roots[index]


def find_roots_bounded(r_new, s_new):
    global r, s
    r = r_new
    s = s_new
    start_guess = [0.5,3]
    a, b = least_squares(equations, start_guess,
                                   bounds = ((0, 1), (1, np.inf)),
                                   gtol=2.220446049250313e-16).x
    return [a,b]


def root_converter(roots):
    a,b = roots
    if a > 1:
        a = 1
    if a < 0:
        a = 0
    if b < 1:
        b = 1

    return [a,b]
