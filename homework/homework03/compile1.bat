@echo off
echo compile problem 1
g++ -c My_Matrix.cpp -O3
g++ -c 1_eigenvalue.cpp -O3 
g++ My_Matrix.o 1_eigenvalue.o -o 1_eigenvalue.exe -O3
del *.o
pause