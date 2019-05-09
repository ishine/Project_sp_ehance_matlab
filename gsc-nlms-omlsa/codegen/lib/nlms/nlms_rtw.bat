@echo off

set MATLAB=D:\matlab

cd .

if "%1"=="" ("D:\matlab\bin\win64\gmake"  -f nlms_rtw.mk all) else ("D:\matlab\bin\win64\gmake"  -f nlms_rtw.mk %1)
@if errorlevel 1 goto error_exit

exit /B 0

:error_exit
echo The make command returned an error of %errorlevel%
An_error_occurred_during_the_call_to_make
