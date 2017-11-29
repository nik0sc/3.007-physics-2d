import numpy as np
# import sympy as sp
import csv
import itertools as it
import functools as ft
from util import *

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
    
def find_total_icm(snapshot, steps, print_lengths=False):
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
    
    if print_lengths:
        for a, b in link_between:
            print('Length between {0} and {1}: {2}'.format(
                a,
                b,
                distance(snapshot[a], snapshot[b])
            ))
    
    snapshot_icm = sum(find_rod_icm(snapshot[a], snapshot[b], steps)
                       for a, b in link_between)
    
    return snapshot_icm

def report(filename, steps, line_dens, break_at=None):
    '''Open filename as csv, calculate moment of inertia with number of steps per link, multiplying the whole thing by the line density line_dens.
    '''
    
    with open(filename) as f:
        coords = read_csv(f)
    icms = []
    print_lengths = True
    for i, snapshot in enumerate(coords):
        if break_at is not None and i == break_at:
            break
        
        icms.append(find_total_icm(snapshot, steps, print_lengths))
        print_lengths = False
        
        if i % 10 == 0:
            print("Finished snapshot index {0}".format(i))
    
    icms = np.array(icms)
    icms *= line_dens
    
    return icms
