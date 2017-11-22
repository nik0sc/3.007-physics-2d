import numpy as np
# import sympy as sp
import csv
import itertools as it
import functools as ft

debug = True

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

def find_rod_icm(start, end, steps):
    '''Find the moment of inertia of a rod *from the origin*, assuming line density = 1.
    
    start: starting coordinate of the rod as tuple
    end: ending coordinate of the rod as tuple
    steps: number of subdivisions of the rod to make
    '''
    # length of line element
    delta_s = distance(start, end) / steps
    
    # divide the line up
    element_genr = divide_line(start, end, steps)
    
    # with each element, sum the square of distances (which is just the sum of squares of coordinate componenents)
    icm = 0
    for element in element_genr:
        icm += np.sum(np.square(element))
    
    icm *= delta_s
    
    assert icm >= 0
    
    return icm
    
def find_total_icm(snapshot, steps):
    '''Find the total moment of inertia for one snapshot in time.
    
    snapshot should be a two-dimensional array of coordinates p_1 to p_7.
    
    LINKS BETWEEN: 
    - origin and 1, 1 and 2, 1 and 6, 2 and 3, 3 and 4, 2 and 4, 3 and 5, 4 and 6, 5 and 6, 5 and 7, 6 and 7
    '''
    link_between = [(0, 1), (1, 2), (1, 6), (2, 3), (3, 4), (2, 4),
                    (3, 5), (4, 6), (5, 6), (5, 7), (6,7)]
    
    assert len(link_between) == 11
    
    # insert origin of the link
    # potential mutability errors
    snapshot = np.insert(snapshot, 0, np.array([[0, 0]]), axis=0)
    
    if debug:
        for a, b in link_between:
            print('Length between {0} and {1}: {2}'.format(
                a,
                b,
                distance(snapshot[a], snapshot[b])
            ))
    
    snapshot_icm = sum(find_rod_icm(snapshot[a], snapshot[b], steps)
                       for a, b in link_between)
    
    return snapshot_icm
                     
