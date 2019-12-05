# USAGE
# Run interface.py to operate with gui

from tkinter import *
from tkinter import messagebox

import rungekutta
import euler
import improved_euler
import exact_solution
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import matplotlib as mpl

"""Combination of matplotlib and tkinter to display the numerical methods for 
solving ODEs. """


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
        elif legline in self.lined_er:
            origline = self.lined_er[legline]
        else:
            origline = self.lined_ger[legline]

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

    def fix_scale(self, vy_ex, vy_er_eu, vy_ger_eu):
        """
        Fix the scale for the graph corresponding to new IVP
        :param vy_ger_eu:
        :param vy_ex: y coordinates of exact solution
        :param vy_er_eu: y coordinates of euler method errors
        :return:
        """
        self.sol_fig.axes[0].set_xlim(left=0, right=self.X)
        self.sol_fig.axes[0].set_ylim(min(vy_ex) - 1, max(vy_ex) + 1)

        self.sol_fig.axes[1].set_xlim(left=0, right=self.X)
        self.sol_fig.axes[1].set_ylim(-0.15 * max(vy_er_eu),
                                      1.5 * max(vy_er_eu))

        self.sol_fig.axes[2].set_xlim(left=self.min_n, right=self.max_n)
        self.sol_fig.axes[2].set_ylim(-1.2 * min(vy_ger_eu),
                                      1.2 * max(vy_ger_eu))

    def initialize_ui_elements(self):
        """
        Initialize various ui elements such as buttons, labels, entries, scale.
        :return: None
        """
        f_top = LabelFrame(text="Vary parameters", bg=self.my_color,
                           font=self.my_font)
        f_top.pack()
        f_bottom = LabelFrame(text="Global Error Interval", bg=self.my_color,
                              font=self.my_font)
        f_bottom.pack()

        padding_label = 20
        padding_entry = 7

        n_text = StringVar()
        n_text.set("Change n value ->")
        n_label = Label(master=f_top, textvariable=n_text, bg=self.my_color,
                        width=padding_label, font=self.my_font)
        n_label.pack(side=LEFT)

        self.n_change = Scale(master=f_top, from_=self.n, to=200,
                              orient=HORIZONTAL, bg=self.my_color,
                              font=self.my_font)
        self.n_change.pack(side=LEFT)

        x0_text = StringVar()
        x0_text.set("Enter x0 value ->")
        x0_label = Label(master=f_top, textvariable=x0_text, bg=self.my_color,
                         width=padding_label, font=self.my_font)

        self.x0_entry = Entry(master=f_top, textvariable=x0_label,
                              bg=self.my_color, font=self.my_font,
                              width=padding_entry)
        self.x0_entry.insert(END, self.x0)

        x0_label.pack(side=LEFT)
        self.x0_entry.pack(side=LEFT)

        y0_text = StringVar()
        y0_text.set("Enter y0 value ->")
        y0_label = Label(master=f_top, textvariable=y0_text, bg=self.my_color,
                         width=padding_label, font=self.my_font)

        self.y0_entry = Entry(master=f_top, textvariable=y0_label,
                              bg=self.my_color, font=self.my_font,
                              width=padding_entry)
        self.y0_entry.insert(END, self.y0)

        y0_label.pack(side=LEFT)
        self.y0_entry.pack(side=LEFT)

        X_text = StringVar()
        X_text.set("Enter X value ->")
        X_label = Label(master=f_top, textvariable=X_text, bg=self.my_color,
                        width=padding_label, font=self.my_font)

        self.X_entry = Entry(master=f_top, textvariable=X_label,
                             bg=self.my_color, font=self.my_font,
                             width=padding_entry)
        self.X_entry.insert(END, self.X)

        X_label.pack(side=LEFT)
        self.X_entry.pack(side=LEFT)

        n_min_text = StringVar()
        n_min_text.set("Min n ->")
        n_min_label = Label(master=f_bottom, textvariable=n_min_text,
                            bg=self.my_color, width=padding_label - 8,
                            font=self.my_font)

        self.n_min_entry = Entry(master=f_bottom, textvariable=n_min_label,
                                 bg=self.my_color, font=self.my_font,
                                 width=padding_entry)
        self.n_min_entry.insert(END, self.min_n)

        n_min_label.pack(side=LEFT)
        self.n_min_entry.pack(side=LEFT)

        n_max_text = StringVar()
        n_max_text.set("Max n ->")
        n_max_label = Label(master=f_bottom, textvariable=n_max_text,
                            bg=self.my_color, width=padding_label - 8,
                            font=self.my_font)

        self.n_max_entry = Entry(master=f_bottom, textvariable=n_max_label,
                                 bg=self.my_color, font=self.my_font,
                                 width=padding_entry)
        self.n_max_entry.insert(END, self.max_n)

        n_max_label.pack(side=LEFT)
        self.n_max_entry.pack(side=LEFT)

        button = Button(master=root, text="Quit", command=self._quit,
                        bg=self.my_color, font=self.my_font)
        button.pack(side=LEFT)
        button = Button(master=root, text="Apply Changes", command=self.update,
                        bg=self.my_color, font=self.my_font)
        button.pack(side=LEFT)

    def __init__(self, root):
        """
        Initialize graphs of the equations, graphs of the local errors.
        Organize initialization of legend. :param root: tkinter root frame
        """
        self.x0 = 1
        self.y0 = 10
        self.X = 15
        self.n = 10
        self.min_n = 10
        self.max_n = 200
        self.my_color = '#F0F0F0'
        self.my_font = ('Roboto Regular', 12)

        plt.style.use('fivethirtyeight')
        mpl.rcParams['lines.linewidth'] = 3
        # initialize solution graphs

        self.runge_kutta = rungekutta.RungeKutta()
        self.euler = euler.Euler()
        self.improved_euler = improved_euler.ImprovedEuler()
        self.exact_solution = exact_solution.ExactSolution()

        vx_rk, vy_rk = self.runge_kutta.get_graph(self.x0, self.y0, self.X,
                                                  self.n)
        vx_eu, vy_eu = self.euler.get_graph(self.x0, self.y0, self.X, self.n)
        vx_ieu, vy_ieu = self.improved_euler.get_graph(self.x0, self.y0, self.X,
                                                       self.n)
        vx_ex, vy_ex = self.exact_solution.exact(self.x0, self.y0, self.X,
                                                 self.n)

        # initialize local error graphs
        vx_er_eu, vy_er_eu = self.euler.local_error(self.x0, self.y0, self.X,
                                                    self.n)
        vx_er_rk, vy_er_rk = self.runge_kutta.local_error(self.x0, self.y0,
                                                          self.X, self.n)
        vx_er_ieu, vy_er_ieu = self.improved_euler.local_error(self.x0, self.y0,
                                                               self.X, self.n)

        vx_ger_rk, vy_ger_rk = self.runge_kutta.total_error(self.x0, self.y0,
                                                            self.X, self.min_n,
                                                            self.max_n)
        vx_ger_eu, vy_ger_eu = self.euler.total_error(self.x0, self.y0, self.X,
                                                      self.min_n, self.max_n)
        vx_ger_ieu, vy_ger_ieu = self.improved_euler.total_error(self.x0,
                                                                 self.y0,
                                                                 self.X,
                                                                 self.min_n,
                                                                 self.max_n)

        # create plots
        self.sol_fig = Figure(figsize=(6, 4), dpi=90)

        sol_subplot = self.sol_fig.add_subplot(221)

        local_error_subplot = self.sol_fig.add_subplot(222)
        local_error_subplot.set_ylabel('Absolute Local Error')
        local_error_subplot.set_xlabel('x')

        total_error_subplot = self.sol_fig.add_subplot(223)
        total_error_subplot.set_ylabel('Total Approximation Error')
        total_error_subplot.set_xlabel('n')

        self.graph_rk, = sol_subplot.plot(vx_rk, vy_rk, label='Runge-Kutta')
        self.graph_eu, = sol_subplot.plot(vx_eu, vy_eu, label='Euler')
        self.graph_ieu, = sol_subplot.plot(vx_ieu, vy_ieu,
                                           label='Improved Euler')
        self.graph_ex, = sol_subplot.plot(vx_ex, vy_ex, label='Exact Solution')
        sol_subplot.set_title("Solution")
        self.leg_sol = sol_subplot.legend(loc='upper left', fancybox=True,
                                          shadow=True)

        self.graph_er_rk, = local_error_subplot.plot(vx_er_rk, vy_er_rk,
                                                     label='Runge-Kutta')
        self.graph_er_eu, = local_error_subplot.plot(vx_er_eu, vy_er_eu,
                                                     label='Euler')
        self.graph_er_ieu, = local_error_subplot.plot(vx_er_ieu, vy_er_ieu,
                                                      label='Improved Euler')
        local_error_subplot.set_title("Local Error")
        self.leg_er = local_error_subplot.legend(loc='upper left',
                                                 fancybox=True, shadow=True)

        self.graph_ger_rk, = total_error_subplot.plot(vx_ger_rk, vy_ger_rk,
                                                      label='Runge-Kutta')
        self.graph_ger_eu, = total_error_subplot.plot(vx_ger_eu, vy_ger_eu,
                                                      label='Euler')
        self.graph_ger_ieu, = total_error_subplot.plot(vx_ger_ieu, vy_ger_ieu,
                                                       label='Improved Euler')
        total_error_subplot.set_title("Total Error")
        self.leg_ger = total_error_subplot.legend(loc='upper left',
                                                  fancybox=True, shadow=True)

        # fix scale corresponding to IVP
        self.fix_scale(vy_ex, vy_er_eu, vy_ger_eu)

        # manage graph picking for solutions
        self.graphs_sol = [self.graph_rk, self.graph_eu, self.graph_ieu,
                           self.graph_ex]
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

        # manage graph picking for global error
        self.graphs_ger = [self.graph_ger_rk, self.graph_ger_eu,
                           self.graph_ger_ieu]
        self.lined_ger = dict()
        for legline, origline in zip(self.leg_ger.get_lines(), self.graphs_ger):
            legline.set_picker(5)  # 5 pts tolerance
            self.lined_ger[legline] = origline

        # put everything on canvas
        self.canvas = FigureCanvasTkAgg(self.sol_fig,
                                        master=root)  # A tk.DrawingArea.
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

        self.toolbar = NavigationToolbar2Tk(self.canvas, root)

        self.toolbar.update()

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
        self.min_n = self.n_min_entry.get()
        self.max_n = self.n_max_entry.get()

        def is_number(s):
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

        # check if IVP is correct
        if x0 != "" and x0 is not None and is_number(x0):
            self.x0 = float(x0)
        else:
            self.show_message("x0 is incorrect")

        if y0 != "" and y0 is not None and is_number(y0):
            self.y0 = float(y0)
        else:
            self.show_message("y0 is incorrect")

        if X != "" and X is not None and is_number(X):
            self.X = float(X)
        else:
            self.show_message("X is incorrect")

        if self.min_n != "" and self.min_n is not None and is_number(
                self.min_n):
            self.min_n = int(self.min_n)
        else:
            self.show_message("Min n is incorrect")

        if self.max_n != "" and self.max_n is not None and is_number(
                self.max_n):
            self.max_n = int(self.max_n)
        else:
            self.show_message("Max n is incorrect")

        if self.min_n >= self.max_n:
            self.show_message(
                "The min n value is larger or equal to the max n value")

        # build graph of solution with new IVP
        vx_rk, vy_rk = self.runge_kutta.get_graph(self.x0, self.y0, self.X, n)
        vx_eu, vy_eu = self.euler.get_graph(self.x0, self.y0, self.X, n)
        vx_ieu, vy_ieu = self.improved_euler.get_graph(self.x0, self.y0, self.X,
                                                       n)
        vx_ex, vy_ex = self.exact_solution.exact(self.x0, self.y0, self.X, n)

        # build graph of local error with new IVP
        vx_er_eu, vy_er_eu = self.euler.local_error(self.x0, self.y0, self.X, n)
        vx_er_rk, vy_er_rk = self.runge_kutta.local_error(self.x0, self.y0,
                                                          self.X, n)
        vx_er_ieu, vy_er_ieu = self.improved_euler.local_error(self.x0, self.y0,
                                                               self.X, n)

        vx_ger_rk, vy_ger_rk = self.runge_kutta.total_error(self.x0, self.y0,
                                                            self.X, self.min_n,
                                                            self.max_n)
        vx_ger_eu, vy_ger_eu = self.euler.total_error(self.x0, self.y0, self.X,
                                                      self.min_n, self.max_n)
        vx_ger_ieu, vy_ger_ieu = self.improved_euler.total_error(self.x0,
                                                                 self.y0,
                                                                 self.X,
                                                                 self.min_n,
                                                                 self.max_n)

        self.fix_scale(vy_ex, vy_er_eu, vy_ger_eu)

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

        self.graph_ger_rk.set_xdata(vx_ger_rk)
        self.graph_ger_rk.set_ydata(vy_ger_rk)

        self.graph_ger_eu.set_xdata(vx_ger_eu)
        self.graph_ger_eu.set_ydata(vy_ger_eu)

        self.graph_ger_ieu.set_xdata(vx_ger_ieu)
        self.graph_ger_ieu.set_ydata(vy_ger_ieu)

        self.canvas.draw()

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
root.geometry("1920x1080")
root.configure(background='#F0F0F0')
app = App(root)
root.mainloop()
