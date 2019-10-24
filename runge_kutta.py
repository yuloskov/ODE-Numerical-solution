import exact_solution
import math


def rk4(f, x0, y0, X, n):
    vx = [0] * (n + 1)
    vy = [0] * (n + 1)
    h = (X - x0) / float(n)
    vx[0] = x = x0
    vy[0] = y = y0
    for i in range(1, n + 1):
        k1 = h * f(x, y)
        k2 = h * f(x + 0.5 * h, y + 0.5 * k1)
        k3 = h * f(x + 0.5 * h, y + 0.5 * k2)
        k4 = h * f(x + h, y + k3)
        vx[i] = x = x0 + i * h
        vy[i] = y = y + (k1 + k2 + k2 + k3 + k3 + k4) / 6
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
        vx[i] = x = x0 + i * h
        k1 = h * f(x_prev, y_prev)
        k2 = h * f(x_prev + 0.5 * h, y_prev + 0.5 * k1)
        k3 = h * f(x_prev + 0.5 * h, y_prev + 0.5 * k2)
        k4 = h * f(x_prev + h, y_prev + k3)
        y = exact_solution.general_solution(c, x)
        vy[i] = math.fabs(y - (y_prev + (k1 + k2 + k2 + k3 + k3 + k4) / 6))
        y_prev = y
        x_prev = x
    return vx, vy