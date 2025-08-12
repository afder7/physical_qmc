import matplotlib.pyplot as plt
import numpy as np

# Создаем фигуру и подграфики 2x2
fig, axs = plt.subplots(2, 4, figsize=(16, 8))
fig.suptitle('8 графика в формате 2×4', fontsize=16)

# Генерируем данные для графиков
x = np.linspace(0, 10, 100)

# График 1: Линейный
axs[0, 0].plot(x, np.sin(x), color='blue')
axs[0, 0].set_title('Синус')
axs[0, 0].grid(True)

# График 2: Квадратичная функция
axs[0, 1].plot(x, x ** 2, color='red')
axs[0, 1].set_title('Квадратичная функция')
axs[0, 1].grid(True)

# График 3: Экспоненциальный рост
axs[0, 2].plot(x, np.exp(x / 5), color='green')
axs[0, 2].set_title('Экспоненциальный рост')
axs[0, 2].grid(True)

# График 4: Случайные данные
random_data = np.random.randn(100).cumsum()
axs[0, 3].plot(x, random_data, color='purple')
axs[0, 3].set_title('Случайные данные')
axs[0, 3].grid(True)


axs[1, 0].plot(x, np.sin(x), color='blue')
axs[1, 0].set_title('Синус')
axs[1, 0].grid(True)

# График 2: Квадратичная функция
axs[1, 1].plot(x, x ** 2, color='red')
axs[1, 1].set_title('Квадратичная функция')
axs[1, 1].grid(True)

# График 3: Экспоненциальный рост
axs[1, 2].plot(x, np.exp(x / 5), color='green')
axs[1, 2].set_title('Экспоненциальный рост')
axs[1, 2].grid(True)

# График 4: Случайные данные
random_data = np.random.randn(100).cumsum()
axs[1, 3].plot(x, random_data, color='purple')
axs[1, 3].set_title('Случайные данные')
axs[1, 3].grid(True)

# Настраиваем отступы между графиками
plt.tight_layout()

# Показываем графики
plt.show()
