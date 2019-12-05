# Usage Get the approximate solution of ODE using Euler's method and
# corresponding errors compare to exact solution

from equation import NumericalSolution


class Euler(NumericalSolution):

    def find_next(self, x0, y0, h):
        return y0 + h * self.f(x0, y0)
