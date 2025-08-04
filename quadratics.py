import matplotlib.pyplot as plt
from random import *


prime = 3
print(prime)

count = 0
x = []
y = []
lin = []
# xr = [randint()]
# yr = []
for i in range(0, prime):
    x.append(i)
    y.append(0)
    for j in range(0, prime):
        if (j ** 2 - i) % prime == 0:
            count += 1
            y[i] = 1
            break
    lin.append(count)


plt.figure(figsize=(12, 6))

plt.scatter(x, y, color="blue")
plt.scatter(x, lin, color="green")
plt.xlabel("Остаток по модулю")
plt.ylabel("Кол-во выбранных чисел <=x")
plt.title(f"Распределение квадратичных вычетов по модулю {prime}")
plt.show()
