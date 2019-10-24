import math


def exact(x0, y0, X, n):
    vx = [0] * (n + 1)
    vy = [0] * (n + 1)
    h = (X - x0) / float(n)
    vx[0] = x0
    vy[0] = y0
    c = solve_ivp(x0, y0)
    for i in range(1, n + 1):
        vx[i] = x = x0 + i * h
        vy[i] = general_solution(c, x)
    return vx, vy


def solve_ivp(x0, y0):
    c = min(math.sqrt(x0) - math.sqrt(y0-x0), math.sqrt(x0) + math.sqrt(y0-x0))
    return c


def general_solution(c, x):
    return 2 * x - 2 * c * math.sqrt(x) + c*c
