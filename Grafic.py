import matplotlib.pyplot as plt
import numpy as np

# для графика в консольку
fig, ax = plt.subplots()
print("Текущие параметры")
l = 60
# int(input("l ="))
T = 60
# int(input("Т ="))
a = 1
# int(input("alfa = "))
k = 0.01
# float(input("k =")
eps = 1e-3
# float(input("epsilon = "))

print('epsilon =', eps)
x = 1
t = 1
#
# def u(x, t):
#     for n in np.array(1, 100):
#         return (16 * (1 - np.cos(np.pi * n)) / (np.pi ** 3 * n ** 3) *
#             np.cos(np.pi * n * a * t / l) * np.sin(np.pi * n * x / l))

def u(x, t):
    res = 0
    for n in range(1, 100):
        res += (16 * (1 - np.cos(np.pi * n)) / (np.pi ** 3 * n ** 3) *
                np.cos(np.pi * n * a * t / l) * np.sin(np.pi * n * x / l))
    return res

xrange = np.arange(0, l, k)
for t in np.arange(0, T + T / 10, 10):
    xs = [u(x, t) for x in xrange]
    plt.ylim(-1, 1)
    plt.xlabel(u't')
    plt.ylabel(u'u(x,t)')
    plt.plot(xrange, xs)
    plt.grid(True)
    plt.savefig('analit_t={}.png'.format(t))
    plt.close()

trange = np.arange(0, T, k)
for x in np.arange(0, l + l / 10, 10):
    xs = [u(x, t) for t in trange]
    plt.ylim(-1, 1)
    plt.xlabel(u'x')
    plt.ylabel(u'u(x,t)')
    plt.plot(trange, xs)
    plt.grid(True)
    plt.savefig('analit_x={}.png'.format(x))
    plt.close()


# Элемент ряда
def u_n(x, t, n):
    return (16 * (1 - np.cos(np.pi * n)) / (np.pi ** 3 * n ** 3) *
            np.cos(np.pi * n * a * t / l) * np.sin(np.pi * n * x / l))


# Сумма ряда для числа элементов elements_num
def u(x, t, elements_num):
    result = 0
    for i in range(1, elements_num + 1):
        result += u_n(x, t, i)
    return result

# Теоретическая оценка количества суммируемых элементов ряда
def N_theoretic(eps):
    return np.ceil((np.sqrt(16 / (np.pi ** 3 * eps)) - 1)).astype(int)

# Эксперементальная оценка количества суммируемых элементов ряда
def redundancy(x, t, eps):
    elements_num = N_theoretic(eps)
    current_u = u(x, t, elements_num)
    for i in range(elements_num, -1, -1):
        if np.abs(current_u - u(x, t, i)) >= eps:
            return i + 1

if __name__ == "__main__":
    print(N_theoretic(eps))
    print(redundancy(1, 1, eps))
    print(u(1, 1, N_theoretic(eps)))
    print(u(1, 1, redundancy(1, 1, eps)))
