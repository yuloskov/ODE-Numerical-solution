# USAGE
# Run interface.py to operate with gui

from tkinter import *
from tkinter import messagebox
import runge_kutta
import euler
import improved_euler
import exact_solution
import equation
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure

"""
Combination of matplotlib and tkinter to display the numerical methods for solving ODEs.
"""


class App:
    def onpick(self, event):
        """
        On the pick event, find the orig line corresponding to the
        legend proxy line, and toggle the visibility

        :param event: user pressed on the graph legend
        :return: None
        """

        legline = event.artist
        if legline in self.lined:
            origline = self.lined[legline]
        else:
            origline = self.lined_er[legline]
        vis = not origline.get_visible()
        origline.set_visible(vis)
        # Change the alpha on the line in the legend so we can see what lines
        # have been toggled
        if vis:
            legline.set_alpha(1.0)
        else:
            legline.set_alpha(0.2)
        self.canvas.draw()

    def on_key_press(self, event):
        """
        Handle events connected with pressing the keys.

        :param event: key pressed
        :return: None
        """

        if (event.key == "enter"):
            self.update()
        key_press_handler(event, self.canvas, self.toolbar)

    def show_message(self, text):
        """
        Show the message to user
        :param text: text to show
        :return: None
        """
        messagebox.showinfo("GUI Python", text)

    def fix_scale(self, vy_ex, vy_er_eu):
        """
        Fix the scale for the graph corresponding to new IVP
        :param vy_ex: y coordinates of exact solution
        :param vy_er_eu: y coordinates of euler method errors
        :return:
        """
        self.sol_fig.axes[0].set_xlim(left=0, right=self.X)
        self.sol_fig.axes[0].set_ylim(min(vy_ex) - 1, max(vy_ex) + 1)

        self.sol_fig.axes[1].set_xlim(left=0, right=self.X)
        self.sol_fig.axes[1].set_ylim(0, 1.5 * max(vy_er_eu))

    def initialize_ui_elements(self):
        """
        Initialize various ui elements such as buttons, labels, entries, scale.
        :return: None
        """
        f_top = LabelFrame(text="Vary parameters")
        f_top.pack()

        self.n_change = Scale(master=f_top, from_=self.n, to=200, label="    Change n:", orient=HORIZONTAL)
        self.n_change.pack(side=BOTTOM)

        x0_text = StringVar()
        x0_text.set("Enter x0 value")
        x0_label = Label(master=f_top, textvariable=x0_text)

        self.x0_entry = Entry(master=f_top, textvariable=x0_label)
        self.x0_entry.insert(END, self.x0)
        self.x0_entry.pack(side=BOTTOM)

        x0_label.pack(side=BOTTOM)

        y0_text = StringVar()
        y0_text.set("Enter y0 value")
        y0_label = Label(master=f_top, textvariable=y0_text)

        self.y0_entry = Entry(master=f_top, textvariable=y0_label)
        self.y0_entry.insert(END, self.y0)
        self.y0_entry.pack(side=BOTTOM)

        y0_label.pack(side=BOTTOM)

        X_text = StringVar()
        X_text.set("Enter X value")
        X_label = Label(master=f_top, textvariable=X_text)

        self.X_entry = Entry(master=f_top, textvariable=X_label)
        self.X_entry.insert(END, self.X)
        self.X_entry.pack(side=BOTTOM)

        X_label.pack(side=BOTTOM)

        button = Button(master=root, text="Quit", command=self._quit)
        button.pack(side=LEFT)
        button = Button(master=root, text="Apply Changes", command=self.update)
        button.pack(side=LEFT)

    def __init__(self, root):
        """
        Initialize graphs of the equations, graphs of the local errors. Organize initialization of legend.
        :param root: tkinter root frame
        """
        self.x0 = 1
        self.y0 = 10
        self.X = 15
        self.n = 10

        # initialize solution graphs
        vx_rk, vy_rk = runge_kutta.rk4(equation.f, self.x0, self.y0, self.X, self.n)
        vx_eu, vy_eu = euler.euler(equation.f, self.x0, self.y0, self.X, self.n)
        vx_ieu, vy_ieu = improved_euler.euler(equation.f, self.x0, self.y0, self.X, self.n)
        vx_ex, vy_ex = exact_solution.exact(self.x0, self.y0, self.X, self.n)

        # initialize local error graphs
        vx_er_eu, vy_er_eu = euler.local_error(equation.f, self.x0, self.y0, self.X, self.n)
        vx_er_rk, vy_er_rk = runge_kutta.local_error(equation.f, self.x0, self.y0, self.X, self.n)
        vx_er_ieu, vy_er_ieu = improved_euler.local_error(equation.f, self.x0, self.y0, self.X, self.n)

        # create plots
        self.sol_fig = Figure(figsize=(5, 4), dpi=100)

        rk_subplot = self.sol_fig.add_subplot(121)
        eu_subplot = self.sol_fig.add_subplot(121)
        ieu_subplot = self.sol_fig.add_subplot(121)
        ex_subplot = self.sol_fig.add_subplot(121)

        er_eu_subplot = self.sol_fig.add_subplot(122)
        er_rk_subplot = self.sol_fig.add_subplot(122)
        er_ieu_subplot = self.sol_fig.add_subplot(122)

        self.graph_rk, = rk_subplot.plot(vx_rk, vy_rk, label='Runge-Kutta')
        self.graph_eu, = eu_subplot.plot(vx_eu, vy_eu, label='Euler')
        self.graph_ieu, = ieu_subplot.plot(vx_ieu, vy_ieu, label='Improved Euler')
        self.graph_ex, = ex_subplot.plot(vx_ex, vy_ex, label='Exact Solution')
        rk_subplot.set_title("Solution")
        self.leg_sol = eu_subplot.legend(loc='upper left', fancybox=True, shadow=True)

        self.graph_er_rk, = er_rk_subplot.plot(vx_er_rk, vy_er_rk, label='Runge-Kutta')
        self.graph_er_eu, = er_eu_subplot.plot(vx_er_eu, vy_er_eu, label='Euler')
        self.graph_er_ieu, = er_ieu_subplot.plot(vx_er_ieu, vy_er_ieu, label='Improved Eurler')
        er_eu_subplot.set_title("Local Error")
        self.leg_er = er_eu_subplot.legend(loc='upper left', fancybox=True, shadow=True)

        # fix scale corresponding to IVP
        self.fix_scale(vy_ex, vy_er_eu)

        # manage graph picking for solutions
        self.graphs_sol = [self.graph_rk, self.graph_eu, self.graph_ieu, self.graph_ex]
        self.lined = dict()
        for legline, origline in zip(self.leg_sol.get_lines(), self.graphs_sol):
            legline.set_picker(5)  # 5 pts tolerance
            self.lined[legline] = origline

        # manage graph picking for errors
        self.graphs_er = [self.graph_er_rk, self.graph_er_eu, self.graph_er_ieu]
        self.lined_er = dict()
        for legline, origline in zip(self.leg_er.get_lines(), self.graphs_er):
            legline.set_picker(5)  # 5 pts tolerance
            self.lined_er[legline] = origline

        # put everything on canvas
        self.canvas = FigureCanvasTkAgg(self.sol_fig, master=root)  # A tk.DrawingArea.
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

        self.toolbar = NavigationToolbar2Tk(self.canvas, root)

        self.toolbar.update()
        self.canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

        self.canvas.mpl_connect("key_press_event", self.on_key_press)
        self.canvas.mpl_connect('pick_event', self.onpick)

        self.initialize_ui_elements()
        mainloop()

    def update(self):
        """
        The function is responsible for updating plots with the change of IVP.
        :return: None
        """
        # get initial values from user
        n = self.n_change.get()
        x0 = self.x0_entry.get()
        y0 = self.y0_entry.get()
        X = self.X_entry.get()

        # check if IVP is correct
        if x0 != "" and x0 is not None and self.is_number(x0):
            self.x0 = float(x0)
        else:
            self.show_message("x0 is incorrect")

        if y0 != "" and y0 is not None and self.is_number(y0):
            self.y0 = float(y0)
        else:
            self.show_message("y0 is incorrect")

        if X != "" and X is not None and self.is_number(X):
            self.X = float(X)
        else:
            self.show_message("X is incorrect")

        # build graph of solution with new IVP
        vx_rk, vy_rk = runge_kutta.rk4(equation.f, self.x0, self.y0, self.X, n)
        vx_eu, vy_eu = euler.euler(equation.f, self.x0, self.y0, self.X, n)
        vx_ieu, vy_ieu = improved_euler.euler(equation.f, self.x0, self.y0, self.X, n)
        vx_ex, vy_ex = exact_solution.exact(self.x0, self.y0, self.X, n)

        # build graph of local error with new IVP
        vx_er_eu, vy_er_eu = euler.local_error(equation.f, self.x0, self.y0, self.X, n)
        vx_er_rk, vy_er_rk = runge_kutta.local_error(equation.f, self.x0, self.y0, self.X, n)
        vx_er_ieu, vy_er_ieu = runge_kutta.local_error(equation.f, self.x0, self.y0, self.X, n)

        self.fix_scale(vy_ex, vy_er_eu)

        # update graphs
        self.graph_rk.set_xdata(vx_rk)
        self.graph_rk.set_ydata(vy_rk)

        self.graph_eu.set_xdata(vx_eu)
        self.graph_eu.set_ydata(vy_eu)

        self.graph_ieu.set_xdata(vx_ieu)
        self.graph_ieu.set_ydata(vy_ieu)

        self.graph_ex.set_xdata(vx_ex)
        self.graph_ex.set_ydata(vy_ex)

        self.graph_er_eu.set_xdata(vx_er_eu)
        self.graph_er_eu.set_ydata(vy_er_eu)

        self.graph_er_rk.set_xdata(vx_er_rk)
        self.graph_er_rk.set_ydata(vy_er_rk)

        self.graph_er_ieu.set_xdata(vx_er_ieu)
        self.graph_er_ieu.set_ydata(vy_er_ieu)

        self.canvas.draw()

    def is_number(self, s):
        """
        Checks if the passed string s is a float number
        :param s: String to be converted to float
        :return: Floating point number converted from s
        """
        try:
            float(s)
            return True
        except ValueError:
            return False

    def _quit(self):
        """
        Quit the app
        :return: None
        """
        root.quit()  # stops mainloop
        root.destroy()  # this is necessary on Windows to prevent
        # Fatal Python Error: PyEval_RestoreThread: NULL tstate


root = Tk()
root.wm_title("Numerical methods for solving ODEs")
app = App(root)
root.mainloop()
