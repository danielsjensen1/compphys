from numpy import fromfunction


def Sfunc(m, n):
    if (m + n) % 2 == 1:
        return 0
    else:
        ans = (2e0 / float(n + m + 5) - 4e0 / float(n + m + 3) + 
               2e0 / float(n + m + 1))
        return ans

def Hfunc(m, n):
    if (m + n) % 2 == 1:
        return 0
    else:
        num = float(1 - m -n - 2 * m * n)
        den = float((m + n + 3) * (m + n + 1) * (m + n - 1))
        return -8e0 * num / den

