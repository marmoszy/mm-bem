# plane wave scattering from soft target in salt water
# input: msh ascii 2.2 file
# output: 1-deg scattering pattern
# MM 26.7.2025            |
from numpy import pi, sqrt, exp, sin, cos, array, zeros
import sys

# input parameters
fname = 'sphere-1.905-600.msh' if len(sys.argv)<=1 else sys.argv[1]
th = 0 if len(sys.argv)<=2 else float(sys.argv[2])
f0 = 38e3 if len(sys.argv)<=3 else float(sys.argv[3])
c0 = 1480 if len(sys.argv)<=4 else float(sys.argv[4])
k = 2*pi*f0/c0; d = [cos(pi/180*th),sin(pi/180*th),0]

# reading nodes and elements from mesh file
with open(fname) as f: lines = f.readlines()
n = int(lines[4])
v = [list(map(float,lines[5+i].split()[1:4])) for i in range(n)]
m = int(lines[n+7])
e = [list(map(int,lines[n+8+i].split()[5:8])) for i in range(m)]

# space - centers of surface triangles
x = [sum([array(v[e[i][j]-1])/3 for j in range(3)]) for i in range(m)]

# Helmholtz single layer operator
S = zeros((m,m),dtype=complex);
for i in range(m):
    for j in range(m):
        r = sqrt(sum((x[i]-x[j])**2))
        S[i][j] = exp(1j*k*r)/r/(4*pi) if r!=0 else 1j*k/(4*pi)

# incident field
f = array([exp(1j*k*p@d) for p in x])

# surface solution
from scipy.linalg import lu_factor, lu_solve
phi = lu_solve(lu_factor(S), -f)

# far field solution
r1 = [[cos(pi/180*i),sin(pi/180*i),0] for i in range(360)]
S = [[exp(-1j*k*p@r)/(4*pi) for p in x] for r in r1]
psc = S @ phi;    
print("\n".join(["%d\t%g"%(i,abs(psc[i])) for i in range(len(psc))]))
