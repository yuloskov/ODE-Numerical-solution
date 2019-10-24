import exact_solution
import math


def euler(f, x0, y0, X, n):
    vx = [0] * (n + 1)
    vy = [0] * (n + 1)
    h = (X - x0) / float(n)
    vx[0] = x = x0
    vy[0] = y = y0
    for i in range(1, n + 1):
        m1 = f(x, y)
        vx[i] = x = x0 + i * h
        m2 = f(x, y + h * m1)
        vy[i] = y = y + h * (m1 + m2) / 2
    return vx, vy


def local_error(f, x0, y0, X, n):
    vx = [0] * (n + 1)
    vy = [0] * (n + 1)
    h = (X - x0) / float(n)
    vx[0] = 0
    vy[0] = 0
    x_prev = x0
    y_prev = y0
    c = exact_solution.solve_ivp(x0, y0)

    for i in range(1, n + 1):
        m1 = f(x_prev, y_prev)
        vx[i] = x = x0 + i * h
        m2 = f(x, y_prev + h * m1)
        y = exact_solution.general_solution(c, x)
        vy[i] = math.fabs(y - (y_prev + h * (m1 + m2) / 2))
        x_prev = x
        y_prev = y
    return vx, vy
