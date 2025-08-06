# MM 19.5.2025: scattering from soft body
import sys
import bempp_cl as bempp
from bempp_cl.api.operators.boundary import helmholtz as boundary
from bempp_cl.api.operators.potential import helmholtz as potential
from bempp_cl.api.linalg import gmres
#from scipy.linalg import gmres
from numpy import pi,sin,cos,exp,abs,array

# input parameters
fname = 'msh/sphere-1.905-600.msh' if len(sys.argv)<=1 else sys.argv[1]
th = 0 if len(sys.argv)<=2 else float(sys.argv[2])
f0 = 38e3 if len(sys.argv)<=3 else float(sys.argv[3])
c0 = 1480 if len(sys.argv)<=4 else float(sys.argv[4])

k0 = 2*pi*f0/c0                             # wave numbers
d = array([cos(pi/180*th),sin(pi/180*th),0])# plane wave direction
#d = array([0, -1.0, 0])                     # plane wave direction

# target mesh and its space
#grid = bempp.api.shapes.sphere(r=0.01905, h=c0/f0/5) # generated mesh
grid = bempp.api.import_grid(fname)
Sp = bempp.api.function_space(grid, "P", 1)

# incident pressure 
@bempp.api.complex_callable
def pinc(x, n, idx, res): res[0] = exp(1j*k0*d.dot(x))
f = bempp.api.GridFunction(Sp, fun=pinc)

# surface solution
S = boundary.single_layer(Sp, Sp, Sp, k0)
u,_ = gmres(S, -f)

# far field solution
#r1 = array([[sin(2*pi*i/360),-cos(2*pi*i/360),0] for i in range(360)]).T
r1 = array([[cos(i*pi/180),sin(i*pi/180),0] for i in range(360)]).T
SL = potential.single_layer(Sp, r1, k0)
psc = - SL.evaluate(u)
print("\n".join(["%d\t%g"%(i,abs(psc[0][i])) for i in range(len(psc[0]))]))


