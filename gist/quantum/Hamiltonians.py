def free_particle(x, p):
    return p * p / 2.


def QHO(x, p):
    return free_particle(x, p) + x * x / 2.


def helium_on_table(x, p):
    V = where(x < 0, -1000, 3) * x
    return free_particle(x, p) + V
