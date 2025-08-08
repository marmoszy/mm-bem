# MM 18.7.2025: scattering from spheres - double transmission problem
import sys
import bempp_cl.api
from bempp_cl.api.operators.boundary import helmholtz, sparse
from bempp_cl.api.operators.potential import helmholtz as potential
from scipy.sparse.linalg import gmres
from numpy import pi,sin,cos,exp,abs,array,block

# input parameters
f0, c0, rho0 = 38e3,    1480, 1024  # sea water
a1, c1, rho1 = 0.01905, 1540, 1045  # fluid bubble
a2, c2, rho2 = 0.01,     340, 1.29  # air bubble

# target mesh and its space
grid = bempp_cl.api.shapes.sphere(r=a1, h=c0/f0/5)  # generated mesh
grid2 = bempp_cl.api.shapes.sphere(r=a2, h=c0/f0/5) # generated mesh
#grid = bempp_cl.api.import_grid("sphere-1.905-600.msh")
#grid2 = bempp_cl.api.import_grid("sphere-1.0-300.msh")
Sp = bempp_cl.api.function_space(grid, "P", 1)
Sp2 = bempp_cl.api.function_space(grid2, "P", 1)
N, N2 = Sp.global_dof_count, Sp2.global_dof_count

# calculated values
k0, k1, k2 = 2*pi*f0/c0, 2*pi*f0/c1, 2*pi*f0/c2 # wave numbers
d = array([1.0, 0, 0])                  # plane wave direction
z = array([0]*N2)

# incident pressure and its normal derivative
@bempp_cl.api.complex_callable
def pinc(x, n, idx, res): res[0] = exp(1j*k0*d.dot(x))
f = bempp_cl.api.GridFunction(Sp, fun=pinc).coefficients
@bempp_cl.api.complex_callable
def dpinc(x, n, idx, res): res[0] = 1j*k0*exp(1j*k0*d.dot(x))*d.dot(n)
g = bempp_cl.api.GridFunction(Sp, fun=dpinc).coefficients

# scaled helmholtz operators
def operators(Sp, Sp2, k, rho):
      S = helmholtz.single_layer(Sp, Sp2, Sp2, k)
      D = helmholtz.double_layer(Sp, Sp2, Sp2, k)
      T = helmholtz.adjoint_double_layer(Sp, Sp2, Sp2, k)
      H = helmholtz.hypersingular(Sp, Sp2, Sp2, k)
      return (rho*rho*S, rho*D, rho*T, -H)
(S0, D0, T0, H0) = operators(Sp, Sp, k0, rho0)
(S1, D1, T1, H1) = operators(Sp, Sp, k1, rho1)
(S2, D2, T2, H2) = operators(Sp2, Sp, k1, rho1)
(S3, D3, T3, H3) = operators(Sp, Sp2, k1, rho1)
(S4, D4, T4, H4) = operators(Sp2, Sp2, k1, rho1)
(S5, D5, T5, H5) = operators(Sp2, Sp2, k2, rho2)
I01 = (rho0+rho1)/2 * sparse.identity(Sp, Sp, Sp)
I12 = (rho1+rho2)/2 * sparse.identity(Sp2, Sp2, Sp2)

# kernel matrix
A = [[ D0 - D1 + I01, S0 - S1,          -D2,           -S2      ],
     [-H0 + H1,      -T0 + T1 + I01,     H2,            T2      ],
     [ D3,                 S3,      D4 - D5 + I12, S4 - S5      ],
     [-H3,                -T3,     -H4 + H5,       T4 - T5 + I12]]
for i in range(len(A)):
      for j in range(len(A[i])):
            A[i][j]=bempp_cl.api.as_matrix(A[i][j].strong_form())

# surface and far field solution
u,_ = gmres(block(A), block([[-f, g/rho0, z, z]]).T)
M = 360               # number of output points on a far field circle
r1 = array([[cos(2*pi*i/M),sin(2*pi*i/M),0] for i in range(M)]).T
SL = potential.single_layer(Sp, r1, k0)
DL = potential.double_layer(Sp, r1, k0)
psi = bempp_cl.api.GridFunction(Sp, coefficients=u[0: N ])
phi = bempp_cl.api.GridFunction(Sp, coefficients=u[N:2*N])
psc = rho0 * DL * psi + rho0**2 * SL * phi                   
print("\n".join(["%d\t%g"%(i,abs(psc[0][i])) for i in range(M)]))
