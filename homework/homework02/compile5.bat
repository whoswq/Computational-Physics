@echo off 
echo compile problem 5
g++ -c integral.cpp -O3
g++ -c 1_D_optimize.cpp -O3
g++ -c 4_berry_phase.cpp -O3
g++ integral.o 1_D_optimize.o 4_berry_phase.o -o 5_berry_phase.exe -O3
del *.o
pause