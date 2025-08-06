# MM 21.5.2025: bempp plane wave scattering from fluid sphere
# python3 fluid-bempp.py fname.msh | gnuplot -p polar.gp
import sys
import bempp_cl as bempp
from bempp_cl.api.operators.boundary import helmholtz as boundary
from bempp_cl.api.operators.potential import helmholtz as potential
from numpy import pi, sin, cos, exp, abs, array, concatenate
from scipy.sparse.linalg import gmres

# input parameters
fname = 'msh/sphere-1.905-600.msh' if len(sys.argv)<=1 else sys.argv[1]
th = 0 if len(sys.argv)<=2 else float(sys.argv[2])
f0 = 38e3 if len(sys.argv)<=3 else float(sys.argv[3])
c0 = 1480 if len(sys.argv)<=4 else float(sys.argv[4])
rho0 = 1024 if len(sys.argv)<=5 else float(sys.argv[5])
c1 = 1540 if len(sys.argv)<=6 else float(sys.argv[6])
rho1 = 1045 if len(sys.argv)<=7 else float(sys.argv[7])

# calculated parameters                      
d = array([cos(pi/180*th),sin(pi/180*th),0])# plane wave direction
k0, k1 = 2*pi*f0/c0, 2*pi*f0/c1             # wave numbers

# target mesh and its space
#grid = bempp.api.shapes.sphere(r=0.01905, h=c0/f0/5) # generated mesh
grid = bempp.api.import_grid(fname)
Sp = bempp.api.function_space(grid, "P", 1)

# incident pressure and its normal derivative
@bempp.api.complex_callable
def pinc(x, n, idx, res): res[0] = exp(1j*k0*d.dot(x))
f = bempp.api.GridFunction(Sp, fun=pinc)
@bempp.api.complex_callable
def dpinc(x, n, idx, res): res[0] = 1j*k0*exp(1j*k0*d.dot(x))*d.dot(n)
g = bempp.api.GridFunction(Sp, fun=dpinc)

# surface solution 
A0 = boundary.multitrace_operator(grid, k0) # external Calderon
A1 = boundary.multitrace_operator(grid, k1) # internal Calderon
A1[0,1] *= rho1/rho0; A1[1,0] *= rho0/rho1  # scaled Calderon
A = (A0+A1).strong_form()
fg = concatenate([f.coefficients, g.coefficients])
u, info = gmres(A*A, A*fg)                  # surface solution

# far field solution
r1 = array([[cos(2*pi*i/360),sin(2*pi*i/360),0] for i in range(360)]).T
SL = potential.single_layer(Sp, r1, k0)
DL = potential.double_layer(Sp, r1, k0)
psi = bempp.api.GridFunction(Sp, coefficients=u[: Sp.global_dof_count])
phi = bempp.api.GridFunction(Sp, coefficients=u[Sp.global_dof_count :])
psc = DL * psi - SL * phi                   # far field solution

# results
print("\n".join(["%d\t%g"%(i,abs(psc[0][i])) for i in range(len(psc[0]))]))
