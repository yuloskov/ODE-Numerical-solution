# Usage
# Get the approximate solution of ODE using Improved Euler's method and corresponding errors compare to exact solution

from exact_solution import ExactSolution
import math
from equation import Equation_Solution


class ImprovedEuler(Equation_Solution):
    f = Equation_Solution.f
    exact_solution = ExactSolution()

    def euler(self, x0, y0, X, n):
        """
        Calculate the approximate solution of given ODE and IVP using Euler's method
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
            m1 = self.f(x, y)
            vx[i] = x = x0 + i * h
            m2 = self.f(x, y + h * m1)
            vy[i] = y = y + h * (m1 + m2) / 2
        return vx, vy

    def local_error(self, x0, y0, X, n):
        """
        Get the local error for Improved Euler's method
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
        c = self.exact_solution.solve_ivp(x0, y0)

        for i in range(1, n + 1):
            m1 = self.f(x_prev, y_prev)
            vx[i] = x = x0 + i * h
            m2 = self.f(x, y_prev + h * m1)
            y = self.exact_solution.general_solution(c, x)
            vy[i] = math.fabs(y - (y_prev + h * (m1 + m2) / 2))
            x_prev = x
            y_prev = y
        return vx, vy

    def global_error(self, x0, y0, X, min_n=10, max_n=100):
        """
        Get the global error of Improved Euler's method for the given equation

        :param f: given equation
        :param x0: initial x
        :param y0: initial y
        :param X: the right side of an interval for which the approximation is calculated
        :param max_n: the maximum number of iterations in the approximation
        :return: vx, vy - the dependence of global error on the number of steps

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
            for i in range(1, cur_n + 1):
                m1 = self.f(x, y)
                x = x0 + i * h
                m2 = self.f(x, y + h * m1)
                y = y + h * (m1 + m2) / 2
            vx.append(cur_n)
            vy.append(math.fabs(self.exact_solution.general_solution(c, x) - y))
        return vx, vy
