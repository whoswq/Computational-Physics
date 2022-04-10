@echo off
echo compile problem 2
g++ -c My_Matrix.cpp -O3
g++ -c 2_eigenvalue.cpp -O3 
g++ My_Matrix.o 2_eigenvalue.o -o 2_eigenvalue.exe -O3
del *.o
pause