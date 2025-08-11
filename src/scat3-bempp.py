# MM 20.7.2025: scattering from three spheres - double transmission problem
import bempp_cl as bempp
from bempp_cl.api.operators.boundary import helmholtz as boundary, sparse
from bempp_cl.api.operators.potential import helmholtz as potential
from scipy.sparse.linalg import gmres
from numpy import pi, sin, cos, exp, abs, array, block
import sys

# --- input parameters
# fname1,fname2,#fname3 = "src/mesh1.msh","src/mesh2.msh","src/mesh3.msh"
fname1 = "msh/sphere-1.905-600.msh"
fname2 = "msh/sphere-0.5-200p.msh"
fname3 = "msh/sphere-0.5-200n.msh"

th0 = 360
f0, c0, rho0 = 38e3,    1480, 1024  # sea water
a1, c1, rho1 = 0.01905, 1540, 1045  # external fluid bubble
a2, c2, rho2 = 0.005,    340, 1.29  # internal air bubble
a3, c3, rho3 = 0.005,    340, 1.29  # internal air bubble

# --- target meshes and its spaces
#grid1 = bempp.api.shapes.sphere(r=a1, h=c0/f0/5, origin=(    0,0,0))
#grid2 = bempp.api.shapes.sphere(r=a2, h=c0/f0/5, origin=( 0.01,0,0)) #  offset
#grid3 = bempp.api.shapes.sphere(r=a3, h=c0/f0/5, origin=(-0.01,0,0)) # -offset
grid1 = bempp.api.import_grid(fname1)
grid2 = bempp.api.import_grid(fname2)
grid3 = bempp.api.import_grid(fname3)

Sp1 = bempp.api.function_space(grid1, "P", 1)
Sp2 = bempp.api.function_space(grid2, "P", 1)
Sp3 = bempp.api.function_space(grid3, "P", 1)

# --- calculated values
k0, k1, k2, k3 = 2*pi*f0/c0, 2*pi*f0/c1, 2*pi*f0/c2, 2*pi*f0/c3
N1,N2,N3 = Sp1.global_dof_count,Sp2.global_dof_count,Sp3.global_dof_count
z2, z3 = array([0]*N2), array([0]*N3)

# --- scaled helmholtz operators
def operators(Sp1, Sp2, k, rho, s=0):
      S = boundary.single_layer(Sp1, Sp2, Sp2, k)
      D = boundary.double_layer(Sp1, Sp2, Sp2, k)
      T = boundary.adjoint_double_layer(Sp1, Sp2, Sp2, k)
      H = boundary.hypersingular(Sp1, Sp2, Sp2, k)
      if s != 0:
            I = s/2 * sparse.identity(Sp1, Sp1, Sp1)
            D, T = D+I, T-I
      return (rho*rho*S, rho*D, -rho*T, H)

print("Assembling ...",file=sys.stderr,flush=True)
(S0, D0, T0, H0) = operators(Sp1, Sp1, k0, rho0,  1)
(S1, D1, T1, H1) = operators(Sp1, Sp1, k1, rho1, -1)
(S2, D2, T2, H2) = operators(Sp2, Sp1, k1, rho1)
(S3, D3, T3, H3) = operators(Sp1, Sp2, k1, rho1)
(S4, D4, T4, H4) = operators(Sp2, Sp2, k1, rho1,  1)
(S5, D5, T5, H5) = operators(Sp2, Sp2, k2, rho2, -1)
(S6, D6, T6, H6) = operators(Sp3, Sp1, k1, rho1)
(S7, D7, T7, H7) = operators(Sp1, Sp3, k1, rho1)
(S8, D8, T8, H8) = operators(Sp3, Sp3, k1, rho1,  1)
(S9, D9, T9, H9) = operators(Sp3, Sp3, k3, rho3, -1)
Z23 = 0 * operators(Sp3, Sp2, k3, rho3)[0]  # N2 x N3 zeros operator
Z32 = 0 * operators(Sp2, Sp3, k3, rho3)[0]  # N3 x N2 zeros operator

# --- kernel matrix
A = [[ D0-D1, S0-S1,  -D2,   -S2,   -D6,    -S6 ],
     [ H0-H1, T0-T1,  -H2,   -T2,   -H6,    -T6 ],
     [   D3,    S3,  D4-D5, S4-S5,  Z23,    Z23 ],
     [   H3,    T3,  H4-H5, T4-T5,  Z23,    Z23 ],
     [   D7,    S7,   Z32,   Z32,  D8-D9, S8-S9 ],
     [   H7,    T7,   Z32,   Z32,  H8-H9, T8-T9 ]]
A = [[bempp.api.as_matrix(x.strong_form()) for x in a] for a in A]

print("Solving ...",file=sys.stderr,flush=True)
for th in range(th0):
      print(th,end=' ',file=sys.stderr,flush=True)
      
      # --- plane wave direction
      d = array([cos(pi/180*th),sin(pi/180*th),0])

      # --- incident pressure and its normal derivative
      @bempp.api.complex_callable
      def pinc(x, n, idx, res): res[0]=exp(1j*k0*d.dot(x))
      f = bempp.api.GridFunction(Sp1, fun=pinc).coefficients
      @bempp.api.complex_callable
      def dpinc(x, n, idx, res): res[0]=1j*k0*exp(1j*k0*d.dot(x))*d.dot(n)
      g = bempp.api.GridFunction(Sp1, fun=dpinc).coefficients

      # --- surface solution
      u,_ = gmres(block(A), block([[-f, g/rho0, z2, z2, z3, z3]]).T)
      psi = bempp.api.GridFunction(Sp1, coefficients=u[0: N1 ])
      phi = bempp.api.GridFunction(Sp1, coefficients=u[N1:2*N1])

      # --- far field solution
      r1 = array([[cos(2*pi*i/360),sin(2*pi*i/360),0] for i in range(360)]).T
      SL = rho0**2 * potential.single_layer(Sp1, r1, k0)
      DL = rho0    * potential.double_layer(Sp1, r1, k0)
      psc = DL * psi + SL * phi                   
      print("\n".join(["%d\t%g"%(i,abs(psc[0][i])) for i in range(len(psc[0]))]))
