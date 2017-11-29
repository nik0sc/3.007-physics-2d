import numpy as np
# import sympy as sp
import csv
import itertools as it
import functools as ft
import util

def snapshot_cm_individual(snapshot):
    '''Find centre of mass of snapshot
    '''
    link_between = [(0, 1), (1, 2), (1, 6), (2, 3), (3, 4), (2, 4),
                (3, 5), (4, 6), (5, 6), (5, 7), (6,7)]
    
    assert len(link_between) == 11
    
    # insert origin of the link
    # potential mutability errors
    snapshot = np.insert(snapshot, 0, np.array([[0, 0]]), axis=0)
    
    # calculate midpoints of links (woo, vectors)
    midpoints = np.array([(snapshot[a] + snapshot[b]) / 2
                         for a, b in link_between])
    
    # calculate "masses" of links (really just distance times a line density
    # factor which is divided away later anyway)
    masses = np.array([util.distance(snapshot[a], snapshot[b]) 
                      for a, b in link_between])
    
    # add it all up by weight...
    assert midpoints.shape == (11, 2), "%s" % repr(midpoints.shape)
    assert masses.shape == (11,), "%s" % repr(masses.shape)
    
    return np.sum(midpoints[i] * masses[i] for i in range(11)) / sum(masses)
