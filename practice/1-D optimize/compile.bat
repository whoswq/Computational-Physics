@echo off
echo first compile library
g++ -c 1_D_optimize.cpp -O3
echo then compile file that need the library
g++ -c test_1_D_optimize.cpp -O3
echo finally link all files
g++ 1_D_optimize.o test_1_D_optimize.o -o test_1_D_optimize.exe -O3
echo delete all .o files
del *.o
pause 