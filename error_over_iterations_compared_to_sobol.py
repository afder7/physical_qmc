import random
from math import *
from scipy.stats import qmc
import matplotlib.pyplot as plt
import copy
from functools import reduce
from operator import mul


number_of_points = 70
dimension = 1000
iter_count = 300
skip_iterations = 10
delta_t = 1 / number_of_points
kappa = 200


eps = 1
m = 1  # actual value is calculated when the points are scattered
c = 1  # actual value is calculated when the points are scattered
A0 = m / delta_t ** 2 + c / (2 * delta_t)
A1 = 2 * m / delta_t ** 2
A2 = m / delta_t ** 2 - c / (2 * delta_t)

cur_u = 0
counter = 0
x_prev = []
x_0 = []

l2_discrepancy_list = []

f = [[]]
it = 0

random_l2_disc = 10000
best_l2_disc = 10000
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
    num = cur_number
    dim = dimension

    U = 0
    for i in range(num):
        for j in range(i + 1, num):
            U += 1 / dist(points[i], points[j])

    return U


def force(points):
    num = cur_number
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
                    re[i][k] += -eps * a_ijk / (d ** 3)

    return re


class PointGenerator:
    def __init__(self, num):
        self.num = num
        self.dim = dimension
        self.points = []
        self.velocities = []

        self.create_points()

    def create_points(self):
        global x_prev, x_0, kappa, m, c, cur_u, l2_discrepancy_list

        for i in range(self.num):
            point = [random.random() for _ in range(self.dim)]
            self.points.append(point)
            x_prev.append(point)
            x_0.append(point)

        global random_l2_disc, best_l2_disc, best_points
        save_points(self.points, "random")
        random_l2_disc = qmc.discrepancy(self.points, method='L2-star')
        l2_discrepancy_list.append(random_l2_disc)
        if random_l2_disc < best_l2_disc:
            best_points = copy.deepcopy(self.points)
            best_l2_disc = random_l2_disc

        u = potential_energy(x_0)

        scalar = 0
        for x_i in x_0:
            for x_ik in x_i:
                scalar += x_ik ** 2
        m = 0.25 * (1 + kappa) * abs(u) * delta_t ** 2 / scalar
        c = -abs(u) * sqrt(kappa) * delta_t / scalar
        cur_u = u

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

    def update(self):
        global x_prev, f

        if it % skip_iterations == 0:
            f = force(self.points)
        else:
            f = [[y / skip_iterations for y in x] for x in f]

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

        if it % skip_iterations == 0:
            global best_l2_disc, best_points, l2_discrepancy_list
            l2_disc = qmc.discrepancy(self.points, method='L2-star')
            l2_discrepancy_list.append(l2_disc)
            if l2_disc < best_l2_disc:
                best_points = copy.deepcopy(self.points)
                best_l2_disc = l2_disc


def save_points(points, save_type):
    with open(f"generation/{save_type}_num_{number_of_points}_dim_{dimension}.txt", "w") as file:
        for pt in points:
            for cd in pt:
                if cd > 1:
                    file.write(str(cd - int(cd)) + " ")
                elif cd < 0:
                    file.write(str(cd - int(cd) + 1) + " ")
                else:
                    file.write(str(cd) + " ")
            file.write("\n")


b_int1 = []
b_int2 = []
b_int3 = []
b_int4 = []

s_int1 = []
s_int2 = []
s_int3 = []
s_int4 = []

b_disc1 = []
s_disc1 = []

for cur_number in range(2, number_of_points + 1):
    app = PointGenerator(cur_number)
    for _ in range(iter_count):
        app.update()
        it += 1

    sobol_engine = qmc.Sobol(d=dimension)  # scramble=True для рандомизации
    sobol_points = list(sobol_engine.random(cur_number))

    d = dimension
    f1 = lambda x: reduce(mul, [1 + sqrt(12) * (x[j] - 0.5) / d for j in range(d)])
    f2 = lambda x: sqrt(12 / d) * (sum(x) - d / 2)
    f3 = lambda x: sqrt(45 / (4 * d)) * (sum([j ** 2 for j in x]) - d / 3)
    f4 = lambda x: sqrt(18 / d) * (sum([sqrt(j) for j in x]) - 2 * d / 3)

    b_int1.append(log(abs(sum([f1(pt) for pt in best_points]) / cur_number), e))
    b_int2.append(log(abs(sum([f2(pt) for pt in best_points]) / cur_number), e))
    b_int3.append(log(abs(sum([f3(pt) for pt in best_points]) / cur_number), e))
    b_int4.append(log(abs(sum([f4(pt) for pt in best_points]) / cur_number), e))

    s_int1.append(log(abs(sum([f1(pt) for pt in sobol_points]) / cur_number), e))
    s_int2.append(log(abs(sum([f2(pt) for pt in sobol_points]) / cur_number), e))
    s_int3.append(log(abs(sum([f3(pt) for pt in sobol_points]) / cur_number), e))
    s_int4.append(log(abs(sum([f4(pt) for pt in sobol_points]) / cur_number), e))

    b_disc1.append(qmc.discrepancy(best_points, method='L2-star'))
    s_disc1.append(qmc.discrepancy(sobol_points, method='L2-star'))

    best_points = []
    best_l2_disc = 10000
    it = 0
    print(cur_number)


fig, axs = plt.subplots(2, 4, figsize=(16, 8))
fig.suptitle('Comparing Sobol and Generated points error over iteration', fontsize=16)

x = list(range(2, number_of_points + 1))

b1 = max(max(s_int1), max(b_int1))
axs[0, 0].plot(x, b_int1, color='blue')
axs[0, 0].set_title("generated 1")
axs[0, 0].set_ylim(-b1, b1)
axs[0, 0].grid(True)

axs[0, 1].plot(x, s_int1, color='blue')
axs[0, 1].set_title("sobol 1")
axs[0, 1].set_ylim(-b1, b1)
axs[0, 1].grid(True)

b2 = min(min(s_int2), min(b_int2))
axs[0, 2].plot(x, b_int2, color='red')
axs[0, 2].set_title("generated 2")
axs[0, 2].set_ylim(b2, 0)
axs[0, 2].grid(True)

axs[0, 3].plot(x, s_int2, color='red')
axs[0, 3].set_title("sobol 2")
axs[0, 3].set_ylim(b2, 0)
axs[0, 3].grid(True)


b3 = min(min(s_int3), min(b_int3))
axs[1, 0].plot(x, b_int3, color='green')
axs[1, 0].set_title("generated 3")
axs[1, 0].set_ylim(b3, 0)
axs[1, 0].grid(True)

# График 2: Квадратичная функция
axs[1, 1].plot(x, s_int3, color='green')
axs[1, 1].set_title("sobol 3")
axs[1, 1].set_ylim(b3, 0)
axs[1, 1].grid(True)


b4 = min(min(s_int4), min(b_int4))
axs[1, 2].plot(x, b_int4, color='purple')
axs[1, 2].set_title("generated 4")
axs[1, 2].set_ylim(b4, 0)
axs[1, 2].grid(True)

axs[1, 3].plot(x, s_int4, color='purple')
axs[1, 3].set_title("sobol 4")
axs[1, 3].set_ylim(b4, 0)
axs[1, 3].grid(True)


fig1, axs1 = plt.subplots(1, 2, figsize=(12, 6))
fig1.suptitle('Comparing Sobol and Generated points discrepancy over iteration', fontsize=16)

x = list(range(2, number_of_points + 1))

b5 = max(max(b_disc1), max(s_disc1))
axs1[0].plot(x, b_disc1, color='blue')
axs1[0].set_title("generated")
axs1[0].set_ylim(0, b5)
axs1[0].grid(True)

axs1[1].plot(x, s_disc1, color='red')
axs1[1].set_title("sobol")
axs1[1].set_ylim(0, b5)
axs1[1].grid(True)


# Настраиваем отступы между графиками
plt.tight_layout()

# Показываем графики
plt.show()

