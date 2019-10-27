# USAGE
# Get coordinates of exact solution of differential equation
import math
from equation import Equation_Solution


class ExactSolution(Equation_Solution):
    def exact(self, x0, y0, X, n):
        """
        The function calculates points of exact solution of IVP of given equation in the file equation.
        :param x0: initial point x0
        :param y0: initial point y0
        :param X: the right side of an interval
        :param n: number of steps
        :return: vx, vy - lists which contain points of the graph
        """
        vx = [0] * (n + 1)
        vy = [0] * (n + 1)
        h = (X - x0) / float(n)
        vx[0] = x0
        vy[0] = y0
        c = self.solve_ivp(x0, y0)
        for i in range(1, n + 1):
            vx[i] = x = x0 + i * h
            vy[i] = self.general_solution(c, x)
        return vx, vy

    def solve_ivp(self, x0, y0):
        """
        Calculate constant for given IVP.
        :param x0: initial point x0
        :param y0: initial point y0
        :return: constant c
        """
        c = min(math.sqrt(x0) - math.sqrt(y0 - x0), math.sqrt(x0) + math.sqrt(y0 - x0))  #x0 - 1 / (-y0 + math.exp(x0))#
        return c

    def general_solution(self, c, x):
        """
        Calculate solution based on constant.
        :param c: constant for the equation
        :param x: the x-coordinate for which you need to get y
        :return: y coordinate based on x coordinate and c
        """
        return 2 * x - 2 * c * math.sqrt(x) + c * c #1/(c-x)+math.exp(x)
