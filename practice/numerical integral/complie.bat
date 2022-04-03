@echo off
:: this is a .bat demo
echo first compile library
g++ -c integral.cpp
echo then compile file that need the library
g++ -c test_integral.cpp
echo finally link all files
g++ integral.o test_integral.o -o test_integral.exe
echo delete all .o files
del *.o
pause 