# USAGE Get the original equation You can put your own differential equation
# here in the return of function f. Be aware that the equation is supposed to
# be in the form y' = f(x, y). You need to put f(x, y) in the function f. If
# you want to plot the exact solution you also need to modify file exact
# solution. There you should an exact solution and express constant of your
# equation for further calculation of IVP.

from math import sqrt
import math
import exact_solution


class NumericalSolution:
    def __init__(self):
        self.exact_solution = exact_solution.ExactSolution()

    def f(self, x, y):
        """
        Get the right side of ODE
        :param x: x-coordinate
        :param y: y-coordinate
        :return: right side of ODE
        """
        return sqrt(y - x) / sqrt(x) + 1

    def find_next(self, x0, y0, h):
        pass

    def get_graph(self, x0, y0, X, n):
        """
        Calculate the approximate solution of given ODE and IVP :param f:
        given equation :param x0: initial x :param y0: initial y :param X:
        the right side of an interval for which the approximation is
        calculated :param n: number of steps :return: vx, vy - coordinates of
        approximate solution
        """
        vx = [0] * (n + 1)
        vy = [0] * (n + 1)
        h = (X - x0) / float(n)
        vx[0] = x = x0
        vy[0] = y = y0
        for i in range(1, n + 1):
            vy[i] = y = self.find_next(x, y, h)
            vx[i] = x = x0 + i * h
        return vx, vy

    def local_error(self, x0, y0, X, n):
        """
        Get the local error for the approximation method
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

        c = self.exact_solution.solve_ivp(x0, y0)
        x_prev = x0
        y_prev = y0
        for i in range(1, n + 1):
            vx[i] = x = x0 + i * h
            y = self.exact_solution.general_solution(c, x)
            vy[i] = math.fabs(y - self.find_next(x_prev, y_prev, h))
            x_prev = x
            y_prev = y
        return vx, vy

    def total_error(self, x0, y0, X, min_n=10, max_n=100):
        """
                Get the total error for approximation method.

                :param f: given equation :param x0: initial x :param y0:
                initial y :param X: the right side of an interval for which
                the approximation is calculated :param max_n: the maximum
                number of iterations in the approximation :return: vx,
                vy - the dependence of global error on the number of steps

                Args:
                    min_n:
                """
        vx = []
        vy = []

        c = self.exact_solution.solve_ivp(x0, y0)
        for cur_n in range(min_n, max_n):
            x = x0
            y = y0
            h = (X - x0) / float(cur_n)
            errors = []
            for i in range(1, cur_n + 1):
                y = self.find_next(x, y, h)
                x = x0 + i * h
                errors.append(
                    math.fabs(self.exact_solution.general_solution(c, x) - y))
            vx.append(cur_n)
            vy.append(max(errors))
        return vx, vy
