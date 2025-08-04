import random
from math import *
from scipy.stats import qmc
import matplotlib.pyplot as plt
import signal
import copy


number_of_points = 4096
dimension = 2
iter_count = 2000
delta_t = 1 / number_of_points
kappa = 200

it = 0

G = 1
m = 1  # actual value is calculated when the points are scattered
c = 1  # actual value is calculated when the points are scattered
A0 = m / delta_t ** 2 + c / (2 * delta_t)
A1 = 2 * m / delta_t ** 2
A2 = m / delta_t ** 2 - c / (2 * delta_t)

cur_U = 0

counter = 0
x_prev = []
x_0 = []

potential_energy_list = []
l2_discrepancy_list = []
wd_discrepancy_list = []
star_discrepancy_list = []

f = [[]]
skip_iterations = 15

best_l2 = 10000
best_points = []


def dist(x_i, x_j):
    dim = len(x_i)
    d = 0
    for k in range(dim):
        d += ((x_i[k] - x_j[k]) * (1 - abs(x_i[k] - x_j[k]))) ** 2
    return sqrt(d)


# points is a set of num points in a given dim
# U_{1,1}
def potential_energy(points):
    num = number_of_points
    dim = dimension

    U = 0
    for i in range(num):
        for j in range(i + 1, num):
            U += 1 / dist(points[i], points[j])

    return U


def force(points):
    num = number_of_points
    dim = dimension
    re = [[0 for __ in range(dim)] for _ in range(num)]

    for i in range(num):
        for j in range(num):
            if j != i:
                d = dist(points[i], points[j])
                for k in range(dim):
                    sign = -1 if points[i][k] < points[j][k] else 1
                    delta = abs(points[i][k] - points[j][k])
                    a_ijk = sign * delta * (1 - delta) * (1 - 2 * delta)
                    re[i][k] += -G * a_ijk / (d ** 3)

    return re


def star_discrepancy(points):
    num = number_of_points
    dim = dimension

    re = 0
    for i in range(num):
        a = points[i]

        sum_ind = 0
        for i in range(num):
            cur_ind = 1
            for j in range(dim):
                cur_ind *= (points[i][j] <= a[j])
            sum_ind += cur_ind

        vol_a = 1
        for j in range(dim):
            vol_a *= a[j]

        re = max(re, sum_ind / num - vol_a)

    return re


class MovingPoints:
    def __init__(self, point_count=number_of_points, dim=dimension):
        self.num = point_count
        self.dim = dim
        self.points = []
        self.velocities = []

        self.create_points()

    def create_points(self):
        global x_prev, x_0, kappa, m, c, cur_U

        for i in range(self.num):
            point = [random.random() for _ in range(self.dim)]
            self.points.append(point)
            x_prev.append(point)
            x_0.append(point)

        U = potential_energy(x_0)

        scalar = 0
        for x_i in x_0:
            for x_ik in x_i:
                scalar += x_ik ** 2
        m = 0.25 * (1 + kappa) * abs(U) * delta_t ** 2 / scalar
        c = -abs(U) * sqrt(kappa) * delta_t / scalar
        cur_U = U

        # calculating x_1
        global f
        f = force(self.points)
        for i in range(self.num):
            for k in range(self.dim):
                self.points[i][k] = self.points[i][k] + 0.5 * (1 / m) * f[i][k] * delta_t ** 2
                if self.points[i][k] > 1:
                    self.points[i][k] = self.points[i][k] - int(self.points[i][k])
                elif self.points[i][k] < 0:
                    self.points[i][k] = self.points[i][k] - int(self.points[i][k]) + 1
        save_points(self, "random")

    def update(self):
        global x_prev, f
        if it % skip_iterations == 0:
            global best_l2, best_points
            f = force(self.points)
            l2_disc = qmc.discrepancy(self.points, method='L2-star')
            if l2_disc < best_l2:
                best_points = copy.deepcopy(self.points)
                best_l2 = l2_disc

        new_x_prev = self.points
        for i in range(self.num):
            for k in range(self.dim):
                x_j = self.points[i][k]
                self.points[i][k] = (-f[i][k] + A1 * x_j - A2 * x_prev[i][k]) / A0
                if self.points[i][k] > 1:
                    self.points[i][k] = self.points[i][k] - int(self.points[i][k])
                elif self.points[i][k] < 0:
                    self.points[i][k] = self.points[i][k] - int(self.points[i][k]) + 1

        x_prev = new_x_prev.copy()

        # global cur_U
        # cur_U = potential_energy(self.points)
        # l2_disc = qmc.discrepancy(self.points, method='L2-star')
        # wd_disc = qmc.discrepancy(self.points, method='WD')
        # star_disc = star_discrepancy(self.points)
        # potential_energy_list.append(cur_U)
        # l2_discrepancy_list.append(l2_disc)
        # wd_discrepancy_list.append(wd_disc)
        # star_discrepancy_list.append(star_disc)


def save_points(app, type):
    with open(f"generation/{type}_num_{app.num}_dim_{app.dim}.txt", "w") as fi:
        for point in app.points:
            for coord in point:
                if coord > 1:
                    fi.write(str(coord - int(coord)) + " ")
                elif coord < 0:
                    fi.write(str(coord - int(coord) + 1) + " ")
                else:
                    fi.write(str(coord) + " ")
            fi.write("\n")


def build_graph():
    plt.figure(figsize=(12, 7))

    plt.subplot(4, 1, 1)  # 4 строки, 1 столбец, позиция 1
    plt.plot(list(range(it)), potential_energy_list, 'g-', label='U')
    plt.legend()

    plt.subplot(4, 1, 2)  # 4 строки, 1 столбец, позиция 2
    plt.plot(list(range(it)), l2_discrepancy_list, 'm--', label='L2')
    plt.legend()

    plt.subplot(4, 1, 3)  # 4 строки, 1 столбец, позиция 3
    plt.plot(list(range(it)), wd_discrepancy_list, 'm--', label='WD')
    plt.legend()

    plt.subplot(4, 1, 4)  # 4 строки, 1 столбец, позиция 4
    plt.plot(list(range(it)), star_discrepancy_list, 'm--', label='Star')
    plt.legend()

    plt.tight_layout()
    plt.show()


app = MovingPoints(dim=dimension, point_count=number_of_points)


def save_on_exit(signum, frame):
    print(f"\nSaving results of current iteration #{it}")
    save_points(app, "modeled")

    print(best_points == app.points)
    with open(f"generation/best_num_{app.num}_dim_{app.dim}.txt", "w") as fi:
        for point in best_points:
            for coord in point:
                if coord > 1:
                    fi.write(str(coord - int(coord)) + " ")
                elif coord < 0:
                    fi.write(str(coord - int(coord) + 1) + " ")
                else:
                    fi.write(str(coord) + " ")
            fi.write("\n")

    # build_graph()
    exit(0)


signal.signal(signal.SIGINT, save_on_exit)
signal.signal(signal.SIGTERM, save_on_exit)


for _ in range(iter_count):
    if it % 10 == 0:
        print("Iteration        ", it)
    app.update()
    it += 1

save_points(app, "modeled")
with open(f"generation/best_num_{app.num}_dim_{app.dim}.txt", "w") as fi:
    for point in best_points:
        for coord in point:
            if coord > 1:
                fi.write(str(coord - int(coord)) + " ")
            elif coord < 0:
                fi.write(str(coord - int(coord) + 1) + " ")
            else:
                fi.write(str(coord) + " ")
        fi.write("\n")
# build_graph()
