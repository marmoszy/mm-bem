# Soft sphere in sea water - MM 4.6.2025 

# ----- functions
P(n,x) = (n>1) ? ((2*n-1)*x*P(n-1,x)-(n-1)*P(n-2,x))/n : (n==1) ? x : 1
j(n,x) = (n>1) ? (2*n-1)/x*j(n-1,x)-j(n-2,x) : (n==1) ? sin(x)/x/x-cos(x)/x : sin(x)/x
y(n,x) = (n>1) ? (2*n-1)/x*y(n-1,x)-y(n-2,x) : (n==1) ? -cos(x)/x/x-sin(x)/x : -cos(x)/x
h(m,x) = j(m,x)+I*y(m,x)
ss(x,t,a) = a/x*sum[m=0:int(x)+5](2*m+1)*j(m,x)/h(m,x)*P(m,cos(t))
# ----- water
a = ARGC>1?ARGV[1]:0.01905; 
f0 = ARGC>2?ARGV[2]:38e3; 
c0 = ARGC>3?ARGV[3]:1480.;
#------ save
set table #sprintf("out/sphere-%g-%g-%g-%g-gp.txt",100*a,0,f0,co)
set xrange [0:360]
set samples 361
unset polar
plot abs(ss(2*pi*f0/c0*a,x*pi/180,a))
unset table
