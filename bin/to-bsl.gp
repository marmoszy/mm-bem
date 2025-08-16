# plot of extracted backscattering length - MM 13.8.2025 

max(x,y)=x>y?x:y;
set polar; set size square; set grid; unset border; set border polar;
unset xtics; unset ytics; set ttics 0,15,360; set samples 361
tsmin=-60;set rrange [tsmin:-20];
set key noenhanced
set key at screen 0.950,0.925
#set key at 135*pi/180,-10

if(ARGC>0) {
  do for [i=1:ARGC] {
    f=ARGV[i]
#    ff=split(f,'.'); 
#    set table ff[1]."-bsl.txt"
    set table $DD
    plot for [i=0:359] f u 1:2 ev ::(i+180)%360:i:(i+180)%360:i w table
    unset table

    plot $DD u ($1*pi/180):(max(tsmin,20*log10($2))) w l t f
  }
}

