#!/bin/bash -v

gcc --version
julia --version
python3 --version
freefem++-mpi
gnuplot --version

gcc src/soft.c -O3 -ffast-math -flto -o bin/soft
time ./bin/soft msh/sphere-1.905-600.msh > out/sphere-1.905-0-38-1480-c.txt
time julia src/soft.jl msh/sphere-1.905-600.msh > out/sphere-1.905-0-38-1480-jl.txt
time python3 src/soft.py msh/sphere-1.905-600.msh > out/sphere-1.905-0-38-1480-py.txt
time freefem++-mpi -v 0 -f src/soft.edp > out/sphere-1.905-0-38-1480-edp.txt
time python3 src/soft-bempp.py msh/sphere-1.905-600.msh > out/sphere-1.905-0-38-1480-bempp.txt
time gnuplot -c src/soft.gp > out/sphere-1.905-0-38-1480-gp.txt

cd out
gnuplot -p -c ../bin/polar.gp sphere-1.905-0-38-1480*.txt
mv polar.svg ../figs/sphere-1.905-0-38-1480.svg
mv polar.pdf ../figs/sphere-1.905-0-38-1480.pdf

gnuplot -p -c ../bin/polar.gp YFT*.txt
mv polar.svg ../figs/YFT-0-38-1480.svg
mv polar.pdf ../figs/YFT-0-38-1480.pdf
cd ..

# ploting analytic results
gnuplot src/soft.gp | gnuplot -p bin/polar.gp
gnuplot src/fluid.gp | gnuplot -p bin/polar.gp
gnuplot src/shell.gp | gnuplot -p bin/polar.gp
