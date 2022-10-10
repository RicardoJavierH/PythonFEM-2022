import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

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

def comparison_plot(f, u, Omega, filename='tmp.pdf'):
  x = sp.Symbol('x')
  f = sp.lambdify([x], f, modules="numpy")
  u = sp.lambdify([x], u, modules="numpy")
  resolution = 401 # no of points in plot
  xcoor = np.linspace(Omega[0], Omega[1], resolution)
  exact = f(xcoor)
  approx = u(xcoor)
  plt.plot(xcoor, approx)
  #plt.hold('on') deprecated
  plt.plot(xcoor, exact)
  plt.legend(['approximation', 'exact'])
  plt.show()
  plt.savefig(filename)

x = sp.Symbol('x')
psi = [1,x]
Omega = [1,2]
f = 10*(x-1)**2-1
u, c =least_squares(f,psi,Omega)
print(u,c)

comparison_plot(f, u, Omega, filename='tmp.pdf')