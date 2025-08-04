from functools import reduce
from operator import mul
from math import *
from scipy.stats import qmc


def test(n, d, f):
    with open(f"generation/modeled_num_{n}_dim_{d}.txt", "r") as fi:
        points = [[float(y) for y in x.split()] for x in fi.read().split("\n")[:-1]]
        s = 0
        for x in points:
            s += f(x)
        print(f"Generated I*      = {s / n}")

    with open(f"generation/random_num_{n}_dim_{d}.txt", "r") as fi:
        points = [[float(y) for y in x.split()] for x in fi.read().split("\n")[:-1]]
        s = 0
        for x in points:
            s += f(x)
        print(f"Random I*         = {s / n}")

    if d < 21000:
        sobol_engine = qmc.Sobol(d=d, scramble=False)  # scramble=True для рандомизации
        points = list(sobol_engine.random(n))
        s = 0
        for x in points:
            s += f(x)
        print(f"Sobol I*          = {s / n}")

        sobol_engine = qmc.Sobol(d=d, scramble=True)  # scramble=True для рандомизации
        points = list(sobol_engine.random(n))
        s = 0
        for x in points:
            s += f(x)
        print(f"Sobol scramble I* = {s / n}")


n1 = 1000
d1 = 5
n2 = 10000
d2 = 5
n3 = 20000
d3 = 10


print("Integral #1")

n = n1
d = d1
f = lambda x: reduce(mul, [1 + sqrt(12) * (x[j] - 0.5) / d for j in range(d)])
print(f"num = {n}, dim = {d}")
test(n, d, f)

n = n2
d = d2
f = lambda x: reduce(mul, [1 + sqrt(12) * (x[j] - 0.5) / d for j in range(d)])
print(f"num = {n}, dim = {d}")
test(n, d, f)

n = n3
d = d3
f = lambda x: reduce(mul, [1 + sqrt(12) * (x[j] - 0.5) / d for j in range(d)])
print(f"num = {n}, dim = {d}")
test(n, d, f)


print("\nIntegral #2")

n = n1
d = d1
f = lambda x: sqrt(12 / d) * (sum(x) - d / 2)
print(f"num = {n}, dim = {d}")
test(n, d, f)

n = n2
d = d2
f = lambda x: sqrt(12 / d) * (sum(x) - d / 2)
print(f"num = {n}, dim = {d}")
test(n, d, f)

n = n3
d = d3
f = lambda x: sqrt(12 / d) * (sum(x) - d / 2)
print(f"num = {n}, dim = {d}")
test(n, d, f)


print("\nIntegral #3")

n = n1
d = d1
f = lambda x: sqrt(45 / (4 * d)) * (sum([j ** 2 for j in x]) - d / 3)
print(f"num = {n}, dim = {d}")
test(n, d, f)

n = n2
d = d2
f = lambda x: sqrt(45 / (4 * d)) * (sum([j ** 2 for j in x]) - d / 3)
print(f"num = {n}, dim = {d}")
test(n, d, f)

n = n3
d = d3
f = lambda x: sqrt(45 / (4 * d)) * (sum([j ** 2 for j in x]) - d / 3)
print(f"num = {n}, dim = {d}")
test(n, d, f)


print("\nIntegral #4")

n = n1
d = d1
f = lambda x: sqrt(18 / d) * (sum([sqrt(j) for j in x]) - 2 * d / 3)
print(f"num = {n}, dim = {d}")
test(n, d, f)

n = n2
d = d2
f = lambda x: sqrt(18 / d) * (sum([sqrt(j) for j in x]) - 2 * d / 3)
print(f"num = {n}, dim = {d}")
test(n, d, f)

n = n3
d = d3
f = lambda x: sqrt(18 / d) * (sum([sqrt(j) for j in x]) - 2 * d / 3)
print(f"num = {n}, dim = {d}")
test(n, d, f)


print("\nIntegral #5")

n = n1
d = d1
f = lambda x: reduce(mul, [x[j] - 0.5 for j in range(d)])
print(f"num = {n}, dim = {d}")
test(n, d, f)

n = n2
d = d2
f = lambda x: reduce(mul, [x[j] - 0.5 for j in range(d)])
print(f"num = {n}, dim = {d}")
test(n, d, f)

n = n3
d = d3
f = lambda x: reduce(mul, [x[j] - 0.5 for j in range(d)])
print(f"num = {n}, dim = {d}")
test(n, d, f)


print("\nIntegral #6")

n = n1
d = d1
f = lambda x: reduce(mul, [sqrt((15 * e ** 15 + 15) / (13 * e ** 15 + 17)) * (e ** (30 * j - 15) - 1) / (e ** (30 * j - 15) + 1) for j in x])
print(f"num = {n}, dim = {d}")
test(n, d, f)

n = n2
d = d2
f = lambda x: reduce(mul, [sqrt((15 * e ** 15 + 15) / (13 * e ** 15 + 17)) * (e ** (30 * j - 15) - 1) / (e ** (30 * j - 15) + 1) for j in x])
print(f"num = {n}, dim = {d}")
test(n, d, f)

n = n3
d = d3
f = lambda x: reduce(mul, [sqrt((15 * e ** 15 + 15) / (13 * e ** 15 + 17)) * (e ** (30 * j - 15) - 1) / (e ** (30 * j - 15) + 1) for j in x])
print(f"num = {n}, dim = {d}")
test(n, d, f)


print("\nIntegral #7")

n = n1
d = d1
f = lambda x: reduce(mul, [(-2.4 * sqrt(7) * (j - 0.5) + 8 * sqrt(7) * (j - 0.5) ** 3) for j in x])
print(f"num = {n}, dim = {d}")
test(n, d, f)

n = n2
d = d2
f = lambda x: reduce(mul, [(-2.4 * sqrt(7) * (j - 0.5) + 8 * sqrt(7) * (j - 0.5) ** 3) for j in x])
print(f"num = {n}, dim = {d}")
test(n, d, f)

n = n3
d = d3
f = lambda x: reduce(mul, [(-2.4 * sqrt(7) * (j - 0.5) + 8 * sqrt(7) * (j - 0.5) ** 3) for j in x])
print(f"num = {n}, dim = {d}")
test(n, d, f)

print("\nIntegral #8")
f = lambda x: reduce(mul, [2 * sqrt(3) * (j - 0.5) for j in x])

n = n1
d = d1
print(f"num = {n}, dim = {d}")
test(n, d, f)

n = n2
d = d2
print(f"num = {n}, dim = {d}")
test(n, d, f)

n = n3
d = d3
print(f"num = {n}, dim = {d}")
test(n, d, f)

print("\nIntegral #9")

n = n1
d = d1
def f(x):
    re = sqrt(2 / (d * (d - 1)))
    sum = 0
    for i in range(d):
        for j in range(i + 1, d):
            ksi = 1 if x[i] < 1 / 6 or x[i] > 4 / 6 else -1
            sum += ksi
    return re * sum
print(f"num = {n}, dim = {d}")
test(n, d, f)

n = n2
d = d2
def f(x):
    re = sqrt(2 / (d * (d - 1)))
    sum = 0
    for i in range(d):
        for j in range(i + 1, d):
            ksi = 1 if x[i] < 1 / 6 or x[i] > 4 / 6 else -1
            sum += ksi
    return re * sum
print(f"num = {n}, dim = {d}")
test(n, d, f)

n = n3
d = d3
def f(x):
    re = sqrt(2 / (d * (d - 1)))
    sum = 0
    for i in range(d):
        for j in range(i + 1, d):
            ksi = 1 if x[i] < 1 / 6 or x[i] > 4 / 6 else -1
            sum += ksi
    return re * sum
print(f"num = {n}, dim = {d}")
test(n, d, f)

