from matplotlib.pyplot import plot, show
from math import exp

def potential(x):
    return (exp(-2*a*x)-e*exp(-a*x))/(pow(a,2))
