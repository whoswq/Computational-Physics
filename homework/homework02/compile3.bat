@echo off
echo compile problem 3
g++ -c integral.cpp -O3
g++ -c 1_D_optimize.cpp -O3
g++ -c 3_trajectory.cpp -O3
g++ integral.o 1_D_optimize.o 3_trajectory.o -o 3_trajectory.exe -O3
del *.o
pause