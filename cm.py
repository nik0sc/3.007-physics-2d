import numpy as np
# import sympy as sp
import csv
import itertools as it
import functools as ft
import util

def find_snapshot_cm(snapshot):
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
    
def report(filename, break_at=None):
    '''Open filename as csv, calculate centre of mass
    '''
    
    with open(filename) as f:
        coords = util.read_csv(f)
    centres_of_mass = []

    for i, snapshot in enumerate(coords):
        if break_at is not None and i == break_at:
            break
        
        centres_of_mass.append(find_snapshot_cm(snapshot))
        
        if i % 10 == 0:
            print("Finished snapshot index {0}".format(i))
    
    return centres_of_mass
