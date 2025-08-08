# shell (two concentric spheres)- MM 7.6.2025 
# ---- functions
P(n,x) = (n>1) ? ((2*n-1)*x*P(n-1,x)-(n-1)*P(n-2,x))/n : (n==1) ? x : 1
j(n,x) = (n>1) ? (2*n-1)/x*j(n-1,x)-j(n-2,x) : (n==1) ? sin(x)/x/x-cos(x)/x : sin(x)/x
jp(n,x) = n/x*j(n,x)-j(n+1,x)
y(n,x) = (n>1) ? (2*n-1)/x*y(n-1,x)-y(n-2,x) : (n==1) ? -cos(x)/x/x-sin(x)/x : -cos(x)/x
yp(n,x) = n/x*y(n,x)-y(n+1,x)
h(m,x) = j(m,x)+I*y(m,x)
hp(m,x) = jp(m,x)+I*yp(m,x)
# ----- Jared NcNew,..,: Sound scattering from two concentric fluid spheres, JASA (2007)
es1(x01,x21,z0,z2,t)=sum[m=0:int(x01)+5]((-I)**m)*(2*m+1)*P(m,cos(t-pi))*h(m,x01/a1)*\
 (z0*jp(m,x21)*j(m,x01)-z2*jp(m,x01)*j(m,x21))/(-z0*jp(m,x21)*h(m,x01)+z2*hp(m,x01)*j(m,x21))
es2(x01,x11,x12,x22, z1,z2, t)=sum[m=0:int(x01)+5]((-I)**m)*(2*m+1)*P(m,cos(t-pi))*h(m,x01/a1)*\
 ((z2*j(m,x12)*jp(m,x22)-z1*jp(m,x12)*j(m,x22))*(z1*hp(m,x11)*j(m,x01)-   h(m,x11)*jp(m,x01))-\
  (z1*hp(m,x12)*j(m,x22)-z2*h(m,x12)*jp(m,x22))*(   j(m,x11)*jp(m,x01)-z1*jp(m,x11)*j(m,x01)))/\
 ((z1*hp(m,x12)*j(m,x22)-z2*h(m,x12)*jp(m,x22))*(   hp(m,x01)*j(m,x11)-z1*h(m,x01)*jp(m,x11))-\
  (z1*jp(m,x12)*j(m,x22)-z2*j(m,x12)*jp(m,x22))*(   hp(m,x01)*h(m,x11)-z1*h(m,x01)*hp(m,x11)))

# ----- params
#f0=38e3;    rho0=1024.; c0=1480.;    
#a1=0.01905; rho1=1045.; c1=1540.;   
#a2=0.01;    rho2=1.29;  c2= 340.;

a1 = ARGC>1?ARGV[1]:0.01905; 
a2 = ARGC>2?ARGV[2]:0.01; 
f0 = ARGC>3?ARGV[3]:38e3; 
c0 = ARGC>4?ARGV[4]:1480.;
rho0 = ARGC>5?ARGV[5]:1024.;
c1 = ARGC>6?ARGV[6]:1540.;
rho1 = ARGC>7?ARGV[7]:1045.;
c2 = ARGC>8?ARGV[8]:340.;
rho2 = ARGC>9?ARGV[9]:1.29;
# -----  
k0=2*pi*f0/c0; k1=2*pi*f0/c1; k2=2*pi*f0/c2;
ka01=k0*a1; ka11=k1*a1; ka12=k1*a2; ka22=k2*a2;
z1=rho0*c0/rho1/c1;     z2=rho0*c0/rho2/c2

#------ write
set table #"shell-gp.txt"
set xrange [0:360]
set samples 361
unset polar
plot abs(es2(ka01,ka11,ka12,ka22, z1,z2, x*pi/180))
unset table
