import exact_solution
import math


def euler(f, x0, y0, X, n):
    vx = [0] * (n + 1)
    vy = [0] * (n + 1)
    h = (X - x0) / float(n)
    vx[0] = x = x0
    vy[0] = y = y0
    for i in range(1, n + 1):
        vx[i] = x = x0 + i * h
        vy[i] = y = y + h * f(x, y)

    return vx, vy


def local_error(f, x0, y0, X, n):
    vx = [0] * (n + 1)
    vy = [0] * (n + 1)
    h = (X - x0) / float(n)
    vx[0] = 0
    vy[0] = 0
    c = exact_solution.solve_ivp(x0, y0)
    y_prev = y0
    x_prev = x0
    for i in range(1, n + 1):
        vx[i] = x = x0 + i * h
        y = exact_solution.general_solution(c, x)
        vy[i] = math.fabs(y - (y_prev + h * f(x_prev, y_prev)))
        y_prev = y
        x_prev = x

    return vx, vy
