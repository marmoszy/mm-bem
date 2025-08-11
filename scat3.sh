#!/bin/bash -v
# scat3 - scattering from three targets (two inside one)
# MM 9.8.2025

# running python+bempp code
time python3 src/scat3-bempp.py > out/scat3-1.905-600-200-200-bempp.txt
cd out
gnuplot -p -e 'oname="../figs/scat3-1.905-600-200-200"' -c ../bin/polar.gp scat3-1.905-600-200-200-bempp.txt
cd ..

# running matlab+gypsilab - macOS specific matlab path!
#/Applications/MATLAB_R2025a.app/bin/matlab -nodisplay -nosplash -nodesktop -r "run('src/#scat3_gypsilab.m');exit;"
#cd out
#gnuplot -p -e 'oname="../figs/scat3-1.905-600-200-200"' -c ../bin/polar.gp scat3-gypsilab-#bsl.txt
#cd ..