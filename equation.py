# USAGE
# Get the original equation
# You can put your own differential equation here in the return of function f.
# Be aware that the equation is supposed to be in the form y' = f(x, y).
# You need to put f(x, y) in the function f.
# If you want to plot the exact solution you also need to modify file exact solution.
# There you should an exact solution and express constant of your equation for further calculation of IVP.

from math import sqrt


class Equation_Solution():
    def f(self, x, y):
        """
        Get the right side of ODE
        :param x: x-coordinate
        :param y: y-coordinate
        :return: right side of ODE
        """
        return sqrt(y - x) / sqrt(x) + 1#math.exp(2*x)+math.exp(x)+y*y-2*y*math.exp(x)#
