@echo off
echo compile latex file
latexmk -xelatex
latexmk -c
pause