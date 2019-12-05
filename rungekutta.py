# Usage Get the approximate solution of ODE using Runge-Kutta method and
# corresponding errors compare to exact solution

from equation import NumericalSolution


class RungeKutta(NumericalSolution):

    def find_next(self, x0, y0, h):
        k1 = h * self.f(x0, y0)
        k2 = h * self.f(x0 + 0.5 * h, y0 + 0.5 * k1)
        k3 = h * self.f(x0 + 0.5 * h, y0 + 0.5 * k2)
        k4 = h * self.f(x0 + h, y0 + k3)
        return y0 + (k1 + k2 + k2 + k3 + k3 + k4) / 6
