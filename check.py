import sympy as sp
import numpy as np
from sympy import nsolve, Tuple
from mpmath import *

# ASSIGNING THEO JANSEN'S 11 HOLY NUMBERS
a = 38.0 * (24/38)
b = 41.5* (24/38)
c = 39.3* (24/38)
d = 30.0* (24/38)
e = 55.8* (24/38)
f = 39.4* (24/38)
g = 30.0* (24/38)
h = 65.7* (24/38)
i = 49.0* (24/38)
j = 50.0* (24/38)
k = 61.9* (24/38)
l = 7.8* (24/38)
m = 15.0* (24/38)

# ASSIGNING COORDINATES
X0 = 0
Y0 = 0
X2 = sp.symbols("X2")
X3 = sp.symbols("X3")
X4 = - 38.0
X5 = sp.symbols("X5")
X6 = sp.symbols("X6")
X7 = sp.symbols("X7")
X8 = 0
Y2 = sp.symbols("Y2")
Y3 = sp.symbols("Y3")
Y4 = - 7.8
Y5 = sp.symbols("Y5")
Y6 = sp.symbols("Y6")
Y7 = sp.symbols("Y7")
Y8 = - 7.8

#DEFINING FUNCTION THAT CALACULATES UNKNOWN POINT FROM 2 KNOWN POINTS; 6 PARAMETERS (REFER TO DIAGRAM IN REPORT)
def solve_coord(XP, YP, L1, XQ, YQ, L2):
    Y = sp.symbols("Y")
    X = sp.symbols("X")
    l = sqrt((YP - YQ)**2 + (XP - XQ)**2)
    m = (YQ - YP)/(XQ - XP)
    tan_alpha = sqrt((2*(L1)*l)**2 - ((L1**2) + l**2 - (L2)**2)**2)/((L1)**2 + l**2 - (L2)**2)
    tan_beta = sqrt((2*(L2)*l)**2 - ((L2**2) + l**2 - (L1)**2)**2)/((L2)**2 + l**2 - (L1)**2)
    tan_tetha = (m +  tan_alpha)/(1- m*(tan_alpha))
    tan_phi = ((m - tan_beta)/(1 + m*tan_beta))
    eqn = Tuple(Y - YP - (tan_tetha)*(X - XP),
                Y - YQ - (tan_phi) * (X - XQ))

    solution = sp.solve(eqn, X, Y, dict = True)
    return solution[0][X].n(), solution[0][Y].n()

#calling coordinates of point by using function as defined above
def coordinates_foot(t):
    X1 = 15 * cos(t)
    Y1 = 15 * sin(t)
    X2, Y2 = solve_coord(X4, Y4, b, X1, Y1, j)
    X3, Y3 = solve_coord(X4, Y4, d, X2, Y2, e)
    X6, Y6 = solve_coord(X1, Y1, k, X4, Y4, c)
    X5, Y5 = solve_coord(X6, Y6, g, X3, Y3, f)
    X7, Y7 = solve_coord(X6, Y6, i, X5, Y5, h)
    return ([X7,Y7])

#iterating result for 1000 timesteps in 1 second
# 0, 0.1, 0.2... 0.9
domain = [t/10 for t in range(0,10)]
result = []

for t in domain:
    X1 = 15 * cos(2 * pi * t)
    Y1 = 15 * sin(2 * pi * t)
    result.append(coordinates_foot(2*pi*t))

print(result)

#c_array = []
#for x in result:
    #complex = x[0] + x[1]*1j
    #c_array.append(complex)


#total_length = sum(abs(p[1]-p[0]) for p in zip(c_array, c_array[1:]))

#print(total_length)




#write coordinates to csv file

# import csv
# fl = open('coordinates_1000_final.csv', 'w')
# writer = csv.writer(fl)
# writer.writerow(['X7', 'Y7']) #if needed
#for values in result:
 #   writer.writerow(values)
#fl.close()


# plot coordinates to graph

#import matplotlib.pyplot as plt
# xs = [x[0] for x in result]
# ys = [x[1] for x in result]
# plt.plot(xs, ys)
# plt.show()

