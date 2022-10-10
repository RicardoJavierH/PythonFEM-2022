import sympy as sp
def least_squares(f, psi, Omega):
  N = len(psi) - 1
  A = sp.zeros(N+1, N+1)
  b = sp.zeros(N+1, 1)
  x = sp.Symbol('x')
  for i in range(N+1):
    for j in range(i, N+1):
      A[i,j] = sp.integrate(psi[i]*psi[j],(x, Omega[0], Omega[1]))
      A[j,i] = A[i,j]
    b[i,0] = sp.integrate(psi[i]*f, (x, Omega[0], Omega[1]))
  c = A.LUsolve(b)
  u = 0
  for i in range(len(psi)):
    u += c[i,0]*psi[i]
  return u, c

x = sp.Symbol('x')
psi = [1,x]
Omega = [0,1]
f = 10*(x-1)**2-1
least_squares(f,psi,Omega)