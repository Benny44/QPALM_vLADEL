set MINICONDA_LIB=C:\Users\u0104126\Downloads\miniconda3\Library\lib
set MINICONDA_INCLUDE=C:\Users\u0104126\Downloads\miniconda3\Library\include

set PATH=%PATH:C:\Program Files\Git\usr\bin;=%

if exist build (cd build && DEL /F/Q/S *.* > NUL && cd .. )

mkdir build
mkdir build\debug
mkdir build\lib

cd build\debug

::cmake ..\.. -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release -DCOVERAGE=ON
cmake ..\.. -DCMAKE_GENERATOR_PLATFORM=x64 -DCOVERAGE=OFF -DUNITTESTS=ON -DNONCONVEX=ON -DINTERFACES=OFF
cmake --build . -v --config Release

:: Run the tests
::\tests\run_all_tests.exe
.\bin\Release\run_all_tests.exe


cd ..\..

if errorlevel 1 exit /b 1