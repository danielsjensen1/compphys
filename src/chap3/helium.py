from chap3.hydrogen import Hydrogen


class Helium(Hydrogen):
    def __init__(self, exponents):
        Hydrogen.__init__(self, exponents)
        self.Q = empty
    
    def coulomb(self, alpha_p, alpha_q):
        return 2 * Hydrogen.coulomb(self, alpha_p, alpha_q)
    
    def fill_arrays(self):
        for index in np.ndindex(self.Q.size()):
            self.S[index] = 
    
    def hartree(self, alpha_p, alpha_q, alpha_r, alpha_s):
        num = 2 * pi**(5 / 2e0)
        den = ((alpha_p + alpha_q) * (alpha_r + alpha_s) *
               sqrt(alpha_p + alpha_q + alpha_r + alpha_s))
        return num / den
    
    def kinetic(self, alpha_p, alpha_q):
        return 2 * Hydrogen.coulomb(self, alpha_p, alpha_q)
    