import matplotlib.pyplot as plt
import numpy as np

# для графика в консольку
fig, ax = plt.subplots()
print("Текущие параметры")
l = 60
# int(input("l ="))
T = 30
# int(input("Т ="))
a = 1
# int(input("alfa = "))
k = 0.01
# float(input("k =")
eps = 1e-3
# float(input("epsilon = "))

print('epsilon =', eps)


h_x = 10
h_t = 10

I = int(l/h_x)
K = int(T/h_t)
j = (a * h_t / h_x)**2

V = np.zeros((K, I))

x = np.arange(0, l, h_x)
t = np.arange(0, T, h_t)


def psi(x):
    return 4 * x / l - (2 * x / l) ** 2


V[0, :] = np.array(list(map(psi, x)))
V[1, :] = np.array(list(map(psi, x)))

# для численных вычислений
for k in np.arange(1, K - 1):
    for i in np.arange(1, I - 1):
        V[k + 1, 0] = 0
        V[k + 1, I] = 0
        V[k + 1, i] = 2 * (1 - j) * V[k, i] + j * (V[k, i + 1] + V[k, i - 1]) - V[k - 1, i]


if __name__ == "__main__":

    mini = min(V)
    maxi = max(V)

    # вывод численных подсчетов
    ch_xrange = np.arange(0, l, I)
    for i in np.arange(0, len(t)):
        ct = [V[i, I] for i in ch_xrange]
        plt.plot(ch_xrange, ct)
        plt.ylim(mini, maxi)
        plt.xlabel(u't')
        plt.ylabel(u'u(x,t)')
        plt.grid(True)
        plt.savefig('test_t={}.png'.format(t))
        plt.close()
    # =========================================

    # trange = np.arange(0, T, k)
    # # -\\-chislenn
    # ch_trange = np.arange(0, T, K)
    #
    # for x in np.arange(0, l + l / 10, 10):
    #     xs = [u(x, t) for t in trange]
    #     plt.ylim(-1, 1)
    #     plt.xlabel(u'x')
    #     plt.ylabel(u'u(x,t)')
    #     # вывод графика, который численно посчитан
    #     plt.plot(x, V(i,))
    #     # -//-
    #     plt.plot(trange, xs)
    #     plt.grid(True)
    #     plt.savefig('analit_x={}.png'.format(x))
    #     plt.close()
    #
    #
    # for i in (1, l / 60):
    #     plt.plot()
    #         x,
    #         V[i,],
    #         type='l',
    #         xlab='x',
    #         ylab='u',
    #         main=paste0('t=', round(t[i], digits=3)),
    #         ylim=c(mini, maxi)
    #
    #
    #     for (i in seq(1, length(x), by=round(length(x) / 60, digits=0))) {
    #     plot(
    #         t,
    #         V[, i],
    #         type = 'l',
    #         xlab = 't',
    #         ylab = 'u',
    #         main = paste0('x=', round(x[i], digits = 3)),
    #         ylim = c(mini, maxi)
    #         )
