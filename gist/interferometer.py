from pylab import *

source = array([0, 0, 5])
wavelngth = 1

detector1 = array([5, 0])
detector2 = array([6, 0])

def get_reading(source, detector):
    return sin(norm(source[:-1] - detector) * 2 * pi / wavelngth)


