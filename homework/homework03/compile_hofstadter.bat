@echo off
echo compile problem Hofstadter
g++ -c My_Matrix.cpp -O3
g++ -c 6_Hofstadter.cpp -O3 
g++ My_Matrix.o 6_Hofstadter.o -o 6_Hofstadter.exe -O3
del *.o
pause