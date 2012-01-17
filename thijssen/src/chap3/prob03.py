from numpy import dot, eye#, sqrt
from numpy.linalg import eigh


def part_a(S):
    """
    Generate a matrix `V` that brings `S` to unit form.
    
    The matrix `V` is the similarity matrix that can be used to create
    a similarity transformation between the overlap matrix `S` and the
    identity matrix `eye(S.size)`.  Symbolically this is written as
    :math:`V^{\dagger} S V = I`.
    
    Parameter(s)
    ------------
    S : array_like, shape (M, M)
         The positive definite Hermitian (or real symmetric) overlap array.
    """
    eigvals, eigvecs = eigh(S)
    U = eigvecs
    #  Create the diagonal matrix containing the inverse square root of
    #  the eigenvalues of S.
    sminushalf = eigvals**(-0.5e0) * eye(eigvals.size)
    V = dot(U, sminushalf)
    return V

def part_b(H, S):
    V = part_a(S)
    #  V dagger (Vdag) is the transpose conjugate of a 
    Vdag = V.conjugate().transpose()
