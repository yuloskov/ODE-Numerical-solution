# Usage
# Get the approximate solution of ODE using Runge-Kutta method and corresponding errors compare to exact solution

import exact_solution
import math


def rk4(f, x0, y0, X, n):
    """
    Calculate the approximate solution of given ODE and IVP using Runge-Kutta method
    :param f: given equation
    :param x0: initial x
    :param y0: initial y
    :param X: the right side of an interval for which the approximation is calculated
    :param n: number of steps
    :return: vx, vy - coordinates of approximate solution
    """
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
    """
    Get the local error for Runge-Kutta method
    :param f: given equation
    :param x0: initial x
    :param y0: initial y
    :param X: the right side of an interval for which the approximation is calculated
    :param n: number of steps
    :return: vx, vy - coordinates of local error
    """
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


def global_error(f, x0, y0, X, min_n=10, max_n=100):
    """
    Get the global error of Runge-Kutta method for the given equation

    :param f: given equation
    :param x0: initial x
    :param y0: initial y
    :param X: the right side of an interval for which the approximation is calculated
    :param max_n: the maximum number of iterations in the approximation
    :return: vx, vy - the dependence of global error on the number of steps
    """
    vx = []
    vy = []

    c = exact_solution.solve_ivp(x0, y0)
    for cur_n in range(min_n, max_n):
        x = x0
        y = y0
        for i in range(1, cur_n):
            h = (X - x0) / float(cur_n)
            k1 = h * f(x, y)
            k2 = h * f(x + 0.5 * h, y + 0.5 * k1)
            k3 = h * f(x + 0.5 * h, y + 0.5 * k2)
            k4 = h * f(x + h, y + k3)
            x = x0 + i * h
            y = y + (k1 + k2 + k2 + k3 + k3 + k4) / 6
        vx.append(cur_n)
        vy.append(math.fabs(exact_solution.general_solution(c, x) - y))
    return vx, vy
