# Fluid sphere - MM 5.6.2025, 7.8.2025

P(n,x)  = (n>1) ? ((2*n-1)*x*P(n-1,x)-(n-1)*P(n-2,x))/n : (n==1) ? x : 1
j(n,x)  = (n>1) ? (2*n-1)/x*j(n-1,x)-j(n-2,x) : (n==1) ? sin(x)/x/x-cos(x)/x : sin(x)/x
jp(n,x) = n/x*j(n,x)-j(n+1,x)
y(n,x)  = (n>1) ? (2*n-1)/x*y(n-1,x)-y(n-2,x) : (n==1) ? -cos(x)/x/x-sin(x)/x : -cos(x)/x
yp(n,x) = n/x*y(n,x)-y(n+1,x)
h(m,x)  = j(m,x)+I*y(m,x)
hp(m,x) = jp(m,x)+I*yp(m,x)
es(x,x1,g,t,a) = a/x*sum[m=0:int(x)+5](2*m+1)*P(m,cos(t))*\
 (g*jp(m,x)*j(m,x1)-(x1/x)*jp(m,x1)*j(m,x))/((x1/x)*jp(m,x1)*h(m,x)-g*hp(m,x)*j(m,x1))

a = ARGC>1?ARGV[1]:0.01905; 
f = ARGC>2?ARGV[2]:38e3; 
c0 = ARGC>3?ARGV[3]:1480.;
rho0 = ARGC>4?ARGV[4]:1024.;
c1 = ARGC>5?ARGV[5]:1540.;
rho1 = ARGC>6?ARGV[6]:1045.;

#------
set table # "fsphere-gp.txt"
set xrange [0:360]
set samples 361
unset polar
plot abs(es(2*pi*f/c0*a,2*pi*f/c1*a,rho1/rho0,x*pi/180,a))
unset table
