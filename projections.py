import matplotlib.pyplot as plt

n = 3125
d = 5

with open(f"generation/modeled_num_{n}_dim_{d}.txt", "r") as f:
    M_points_2d = [x.split()[:2] for x in f.read().split("\n")][:-1]
    x1 = [float(i[0]) for i in M_points_2d]
    y1 = [float(i[1]) for i in M_points_2d]

with open(f"generation/random_num_{n}_dim_{d}.txt", "r") as f:
    R_points_2d = [x.split()[:2] for x in f.read().split("\n")][:-1]
    x2 = [float(i[0]) for i in R_points_2d]
    y2 = [float(i[1]) for i in R_points_2d]

plt.figure(figsize=(7, 14))
print(x1)
plt.subplot(2, 1, 1)  # 3 строки, 1 столбец, позиция 1
plt.scatter(x1, y1, label='Generated')
plt.legend()

plt.subplot(2, 1, 2)  # 3 строки, 1 столбец, позиция 2
plt.scatter(x2, y2, label='Random')
plt.legend()

plt.tight_layout()
plt.show()
