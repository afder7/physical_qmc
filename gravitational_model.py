import tkinter as tk
import random
from math import *
import numpy as np
from scipy.stats import qmc

G = 1
delta_t = 0.001
kappa = 100
m = 1  # actual value is calculated when the points are scattered
c = 1  # actual value is calculated when the points are scattered
A0 = m / delta_t ** 2 + c / (2 * delta_t)
A1 = 2 * m / delta_t ** 2
A2 = m / delta_t ** 2 - c / (2 * delta_t)

cur_U = 0

counter = 0
x_prev = []
x_0 = []


class Point:
    def __init__(self, img, coords):
        self.img = img
        self.coords = coords


def dist(x_i, x_j):
    dim = len(x_i)
    d = 0
    for k in range(dim):
        d += ((x_i[k] - x_j[k]) * (1 - abs(x_i[k] - x_j[k]))) ** 2
    return sqrt(d)


# points is a set of num points of class Point in a given dim
# U_{1,1}
def potential_energy(points):
    num = len(points)
    dim = len(points[0].coords)

    U = 0
    for i in range(num):
        for j in range(i + 1, num):
            U -= G * m / dist(points[i].coords, points[j].coords)

    print(U)
    return U


def force(points):
    num = len(points)
    dim = len(points[0].coords)
    f = [[0 for __ in range(dim)] for _ in range(num)]

    for i in range(num):
        for j in range(num):
            if j != i:
                for k in range(dim):
                    sign = -1 if points[i].coords[k] < points[j].coords[k] else 1
                    delta = abs(points[i].coords[k] - points[j].coords[k])
                    a_ijk = sign * delta * (1 - delta) * (1 - 2 * delta)
                    f[i][k] += -G * a_ijk / dist(points[i].coords, points[j].coords) ** 3

    return f


class MovingPoints:
    def __init__(self, root, width=700, height=400, point_count=100, dim=2):
        self.root = root
        self.width = width
        self.height = height
        self.num = point_count
        self.dim = dim
        self.points = []
        self.velocities = []

        self.canvas = tk.Canvas(root, width=width, height=height, bg="white")
        self.canvas.pack()

        self.draw_grid()
        self.create_points()

        self.update()

    def draw_grid(self, spacing=50):
        for x in range(0, self.width, spacing):
            self.canvas.create_line(x, 0, x, self.height, fill="lightgray")
        for y in range(0, self.height, spacing):
            self.canvas.create_line(0, y, self.width, y, fill="lightgray")

    def create_points(self):
        global x_prev, x_0, kappa, m, c, cur_U

        for i in range(self.num):
            x = random.randint(0, self.width)
            y = random.randint(0, self.height)

            img = self.canvas.create_oval(x - 2, y - 2, x + 2, y + 2, fill="blue")
            point = Point(img, coords=[x / self.width, y / self.height])
            self.points.append(point)
            x_prev.append(point)
            x_0.append(point)

        U = potential_energy(x_0)

        scalar = 0
        for x_i in x_0:
            for x_ik in x_i.coords:
                scalar += x_ik ** 2
        m = 0.25 * (1 + kappa) * abs(U) * delta_t ** 2 / scalar
        c = abs(U) * sqrt(kappa) * delta_t / scalar
        cur_U = U

        # calculating x_1
        f = force(self.points)
        for i in range(self.num):
            for k in range(self.dim):
                self.points[i].coords[k] = self.points[i].coords[k] + 0.5 * (1 / m) * f[i][k] * delta_t ** 2

    def update(self):
        global counter, x_prev
        if counter < 1000:

            f = force(self.points)

            new_x_prev = self.points
            for i in range(self.num):
                self.canvas.delete(self.points[i].img)

                for k in range(self.dim):
                    x_j = self.points[i].coords[k]
                    self.points[i].coords[k] = (A1 * x_j - A2 * x_prev[i].coords[k] - f[i][k]) / A0
                    if self.points[i].coords[k] > 1:
                        self.points[i].coords[k] = self.points[i].coords[k] - int(self.points[i].coords[k])
                    elif self.points[i].coords[k] < 0:
                        self.points[i].coords[k] = self.points[i].coords[k] - int(self.points[i].coords[k]) + 1

                x = self.points[i].coords[0] * self.width
                y = self.points[i].coords[1] * self.height
                self.points[i].img = self.canvas.create_oval(x - 2, y - 2, x + 2, y + 2, fill="blue")

            x_prev = new_x_prev.copy()

            self.root.after(int(100000 * delta_t), self.update)
            counter += 1

            global cur_U

            cur_U = potential_energy(self.points)
            point_array = []
            for point in self.points:
                point_array.append(point.coords)
            point_array = np.array(point_array)
            l2_disc = qmc.discrepancy(point_array, method='L2-star')

            print(f"Potential Energy", cur_U)
            print(f"L2 discrepancy  ", l2_disc)


# Create and run the animation
root = tk.Tk()
app = MovingPoints(root)
root.mainloop()