# Usage Get the approximate solution of ODE using Improved Euler's method and
# corresponding errors compare to exact solution

from equation import NumericalSolution


class ImprovedEuler(NumericalSolution):

    def find_next(self, x0, y0, h):
        m1 = self.f(x0, y0)
        x0 += h
        m2 = self.f(x0, y0 + h * m1)
        return y0 + h * (m1 + m2) / 2
