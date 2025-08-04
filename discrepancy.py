from scipy.stats import qmc
from math import *

num = 4096
dim = 2


with open(f"generation/random_num_{num}_dim_{dim}.txt", "r") as f:
    points = []
    k = 0
    for i in f.readlines():
        points.append(list(map(float, i.strip().split())))
        k += 1
    r_l2_disc = qmc.discrepancy(points, method='L2-star')

with open(f"generation/modeled_num_{num}_dim_{dim}.txt", "r") as f:
    points = []
    for i in f.readlines():
        points.append(list(map(float, i.strip().split())))
    m_l2_disc = qmc.discrepancy(points, method='L2-star')

with open(f"generation/best_num_{num}_dim_{dim}.txt", "r") as f:
    points = []
    for i in f.readlines():
        points.append(list(map(float, i.strip().split())))
    b_l2_disc = qmc.discrepancy(points, method='L2-star')


print("Random L2-star discrepancy: ", r_l2_disc)
print("Modeled L2-star discrepancy:", m_l2_disc)
print("Best L2-star discrepancy:   ", b_l2_disc)
print("")

