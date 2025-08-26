f = "mm-bem.txt"

fx(x) = 20*log10(x)
xx(u, v, r) = fx(r) * cos(v*pi/180) * cos(u*pi/180)
yy(u, v, r) = fx(r) * cos(v*pi/180) * sin(u*pi/180)
zz(u, v, r) = fx(r) * sin(v*pi/180)

set view equal xyz
set cbrange [-60:-20]

#set palette defined (0 "black", 0 "black", 0 "dark-blue", 1 "yellow")
#set pm3d
#set pm3d corners2color c1
#set pm3d scansauto border lc "black" lw 0.5
#set contourfill ztics

#splot f u (cos($1*pi/180)*cos($2*pi/180)*20*log10($3)):(cos($1*pi/180)*sin($2*pi/180)*20*log10($3)):(sin($1*pi/180)*20*log10($3)) w pm3d t f 

unset border
unset xtics
unset ytics
unset ztics

splot f u (xx($1,$2,$3)):(yy($1,$2,$3)):(zz($1,$2,$3)):(fx($3)) w pm3d t f, f u (xx($1,$2,$3)):(yy($1,$2,$3)):(zz($1,$2,$3)):(fx($3)) w l lc "black" t ""
