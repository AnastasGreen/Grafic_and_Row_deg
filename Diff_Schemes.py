import matplotlib.pyplot as plt
import numpy as np

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

h_x = 0.5
h_t = 0.5

I = int(l / h_x)
K = int(T / h_t)
j = (a * h_t / h_x) ** 2


V = np.zeros((K, I))

x = np.linspace(0, l, num=I)
t = np.linspace(0, T, num=K)
print(round(len(t)/2))
print(len(x))

def u(x, t):
    res = 0
    for n in range(1, 100):
        res += (16 * (1 - np.cos(np.pi * n)) / (np.pi ** 3 * n ** 3) *
                np.cos(np.pi * n * a * t / l) * np.sin(np.pi * n * x / l))
    return res


def psi(x):
    return 4 * x / l - (2 * x / l) ** 2


V[0, :] = np.array(list(map(psi, x)))
V[1, :] = np.array(list(map(psi, x)))

# для численных вычислений
for k in range(1, K - 1):
    V[k + 1, 0] = 0
    V[k + 1, I - 1] = 0
    for i in np.arange(1, I - 2):
        V[k + 1, i] = 2 * (1 - j) * V[k, i] + j * (V[k, i + 1] + V[k, i - 1]) - V[k - 1, i]

maxi = np.max(V)
mini = np.min(V)

if __name__ == "__main__":

    # вывод численных подсчетов
    # правильный варик
    ch_xrange = np.arange(0, l, h_x)
    for i in range(0, round(len(t)), round(len(t) / 15)):
        ct = V[i, :]
        plt.ylim(mini, maxi)
        plt.xlabel(u't')
        plt.ylabel(u'u(x,t)')
        plt.plot(ch_xrange, ct, label='chislen')
        plt.grid(True)
        plt.savefig('test_t={}.png'.format(t[i]))
        plt.close()
        # конец
    xrange = np.arange(0, l, k)
    for t in np.arange(0, T + T / 10, 10):
        xs = [u(x, t) for x in xrange]
        plt.ylim(-1, 1)
        # plt.xlabel(u't')
        # plt.ylabel(u'u(x,t)')
        plt.plot(xrange, xs, label='analit')
        plt.grid(True)
        plt.savefig('analit_t={}.png'.format(t))
        plt.close()
        # =========================================
        # ct = V[I-1, :]
        # plt.plot(ch_xrange, ct)
        # plt.ylim(mini, maxi)
        # plt.xlabel(u't')
        # plt.ylabel(u'u(x,t)')
        # plt.grid(True)
        # plt.savefig('t={}.png'.format(t[I-1]))
        # plt.close()
    # =========================================

    # правильный вариант

    ch_trange = np.arange(0, T, h_t)
    for i in range(0, round(len(x)), round(len(x) / 15)):

        cx = V[:, i]
        plt.ylim(mini, maxi)
        plt.xlabel(u'x')
        plt.ylabel(u'u(x,t)')
        plt.plot(ch_trange, cx)
        plt.grid(True)
        plt.savefig('test_x={}.png'.format(x[i]))
        plt.close()

    # =========================================

    # ch_trange = np.arange(0, T, h_t)
    # for i in range(0, len(t), round(len(t) / 10)):
    #     cx = V[:, i]
    #     plt.plot(ch_trange, cx)
    #     plt.ylim(mini, maxi)
    #     plt.xlabel(u'x')
    #     plt.ylabel(u'u(x,t)')
    #     plt.grid(True)
    #     plt.savefig('test_x={}.png'.format(x[i]))
    #     plt.close()
    # =========================================

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
