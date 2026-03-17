"""
This code carries out the calculations in the proof of Theorem (thm:Z-ideal-sl3)
and the subsequent Remark, from "Singularities of Steinberg deformation rings" by 
Daniel Funck and Jack Shotton (https://arxiv.org/abs/2409.17812). Code produced with 
the aid of Claude Code.

We work in bidegree (2,1) of S = Sym((sl_3^*)^2), with M, N generic traceless
3x3 matrices over ZZ.  We show:
  (1) the entries of M^2*N and N*M^2 span a 17-dim space mod (S_+ * I_tilde);
  (2) entries of M*N*M lie in I_tilde;
  (3) any two of {M^2*N, N*M^2, M*N*M} give the same 17-dim quotient.
Smith normal forms over ZZ show the only bad prime is 2.
"""

R = PolynomialRing(ZZ, 18, 'x')

M, N = matrix(R, 3, 3, R.gens()[:9]), matrix(R, 3, 3, R.gens()[9:])

# --- Bidegree-(2,1) polynomials ---
M2N, NM2, MNM = M*M*N, N*M*M, M*N*M
N_tr2M = (M.trace())^2 * N
N_trM2 = (M*M).trace() * N   
M_trMN = (M*N).trace() * M  
M2_trN = N.trace() * M^2
MN_trM = M.trace() * M * N
NM_trM = M.trace() * N*M
trM_trN_M = M.trace() * N.trace() * M

# List of spanning set of (2,1)-part of S_+\tilde{I}
lowerpolylist = (N_trM2, M_trMN, N_tr2M, M2_trN, MN_trM, NM_trM, trM_trN_M)
lowerpolys = []
for f in lowerpolylist:
    lowerpolys += f.list()

# Spanning set of \tilde{I} in bidegree (2,1)
all_polys = M2N.list() + NM2.list() + MNM.list() + lowerpolys
# All monomials that occur in these polynomials
monoms = sorted({m for f in all_polys for m in f.monomials()}, key=str)

def vec(f):
    return [f.monomial_coefficient(m) for m in monoms]

def to_mat(polys):
    return matrix(ZZ, [vec(f) for f in polys])

# Column span of W is S_+ \tilde{I} in bidegree (2,1)
W = to_mat(lowerpolys)
# Column span of A is J = S_+ \tilde{I} + (coeffs M^2N, NM^2) in bidegree (2,1)
A = W.stack(to_mat(M2N.list() + NM2.list()))       
# Column span of B is \tilde{I} in bidegree (2,1)
B = A.stack(to_mat(MNM.list()))                 


rk_W, rk_A, rk_B = W.rank(), A.rank(), B.rank()

# Check that bidegree (2,1) part of J/S^+\tilde{I} has dimension 17 (in characteristic 0)
print(f"dim(J/S^+\\tilde{I}) = {rk_A - rk_W}  [expected 17]")
# Check that J agrees with \tilde{I} (in characteristic 0)
print(f"MNM in span: {rk_B == rk_A}  [expected True]")

# --- Smith normal forms: find bad primes ---
def bad_primes(mat):
    eds = mat.elementary_divisors()
    if eds[0] ==0:
       return []
    for i in range(len(eds)-1):
       if eds[i+1] == 0:
       	  return [p for p, _ in factor(eds[i])]
    
# Since 2 is only bad prime for any of W, A, B, above calculation is valid in characteristic > 2
print(f"Bad primes for S_+\\tilde{I}: {bad_primes(W)}")
print(f"Bad primes for J: {bad_primes(A)}")
print(f"Bad primes for \\tilde{I}: {bad_primes(B)}")

# Proof of Remark: any two of {M^2N, NM^2, MNM} suffice
C = W.stack(to_mat(M2N.list() + MNM.list()))
rk_C = C.rank()
print(f"dim(M2N, MNM mod W) = {rk_C - rk_W}  [expected 17]")

print(f"Bad primes: {bad_primes(C)} [expected [2]]")