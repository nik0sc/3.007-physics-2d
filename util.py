import numpy as np
# import sympy as sp
import csv
import itertools as it
import functools as ft

def read_csv(f):
    '''Read coordinate data from file object f and create numpy arrays.
    
    WARNING: Performs eval() on each cell of the csv file
    '''
    
    csv_reader = csv.reader(f)
    # each row is a "snapshot" in time of the linkage state
    # each cell is actually a tuple of coordinates
    arr = [[eval(cell) for cell in row] for row in csv_reader]

    return np.array(arr)


def divide_line(a, b, steps):
    '''Make a generator for midpoints along the line (xa, ya), (xb, yb) in however many steps'''
    # no relation to calculus
    d = b - a
    
    for t in range(steps):
        # t in [0,1]
        normalised_t = t / steps
        assert 0 <= normalised_t < 1
        
        # m = a + t(b-a)
        midpoint = a + d * normalised_t
        
        yield midpoint
    
    
def distance(a, b=None):
    if b is not None:
        return np.linalg.norm(a - b)
    else:
        return np.linalg.norm(a)
