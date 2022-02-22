# optimizeCostFunctionParameter.py>

from runInverseSolver import *
from scipy.optimize import minimize_scalar 
from scipy import optimize

mesh = "20 5 15"
Deltat = 0.5
cleanCase = 0

f = lambda x: runInverseSolver(mesh, Deltat, x, cleanCase)

#res = minimize_scalar(f, bracket=(1e-11, 1e-8),bounds=(1e-16, 1e-6), method='brent')
res = optimize.minimize(f, 1e-6, bounds=[(1e-16, 1e-6)],  method="Nelder-Mead")

print('Optimal Input x: %.2e' % res["x"])
print('Optimal Output f(x): %.2e' % res["fun"])
print('Total Evaluations n: %d' % res['nfev'])
