from numpy import cos, ndenumerate, pi, sin, zeros
from chap3.prob03 import part_b


class SquareWell(object):
    def __init__(self, N=5):
        self.N = N
        S = zeros((N, N))
        H = zeros((N, N))
        for (m, n), val in ndenumerate(S):
            if (m + n) % 2 == 0:
                S[m, n] = (2e0 / float(n + m + 5) - 
                           4e0 / float(n + m + 3) + 
                           2e0 / float(n + m + 1))
                num = float(1 - m -n - 2 * m * n)
                den = float((m + n + 3) * (m + n + 1) * (m + n - 1))
                H[m, n] = -8e0 * num / den
        self.S = S
        self.H = H
    
    def approx_eigfunc(self, n, x):
        C = self.eigvecs[:, n]
        f = 0
        for i in range(0, self.N):
            f += C[i] * self.psi(i, x)
        return f
    
    def exact_eigfunc(self, n, x):
        kn = n * pi / 2e0
        if n % 2 == 1:
            return cos(kn * x)
        else:
            return sin(kn * x)
    
    def psi(self, n, x):
        return x**n * (x - 1) * (x + 1)
    
    def variational(self):
        self.eigvals, self.eigvecs = part_b(self.H, self.S)
