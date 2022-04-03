@echo off
echo compile problem 2
g++ -c integral.cpp -O3
g++ -c 1_D_optimize.cpp -O3
g++ -c 2_trajectory.cpp -O3
g++ integral.o 1_D_optimize.o 2_trajectory.o -o 2_trajectory.exe -O3
del *.o
pause