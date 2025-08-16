# txt log polar multiplot - MM 7.6.2025 
# ex: gnuplot -p -c polar.gp *.txt
# ex: gnuplot -e 'oname="x"' -c polar.gp soft-c.txt soft-py.txt
# output file name: gnuplot -p -e 'oname="all"' -c polar.gp *.txt
# only record 270: gnuplot -p -e "bno=270" -c polar.gp ../out/simplify50_MACK_flesh_L20.txt

#fnames="soft-c.txt soft-py.txt soft-jl.txt soft-js.txt"
fnames='<cat'
if(ARGC>0) {
  fnames=''
  do for [i=1:ARGC] {
    fnames=fnames.' '.ARGV[i]
  }
}
if(!exists("oname")) oname="polar" # output filename

max(x,y)=x>y?x:y;
set polar; set size square; set grid; unset border; set border polar;
unset xtics; unset ytics; set ttics 0,15,360; set samples 361
tsmin=-60;set rrange [tsmin:-20];
set key noenhanced
set key at screen 0.950,0.925
#set key at 135*pi/180,-10

if(exists("bno")) {
plot for [f in fnames] sprintf("%s",f) u ($1*pi/180):(max(tsmin,20*log10($2))) every :::bno::bno w l t f
} else {
plot for [f in fnames] sprintf("%s",f) u ($1*pi/180):(max(tsmin,20*log10($2))) w l t f
}

if(fnames ne '<cat') {
  #replot
  #set term pdfcairo; set output sprintf("%s.pdf",oname);
  #replot
  set term svg; set output sprintf("%s.svg",oname);
  replot
}
