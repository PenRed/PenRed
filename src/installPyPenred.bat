mkdir build
cd build


cmake -DWITH_DICOM="OFF" -DWITH_MULTI_THREADING="ON" -DWITH_MPI="OFF" -DWITH_LB="OFF" -DBUILD_PYTHON_MODULES="ON" -DBUILD_TESTS="OFF" -DBUILD_UTILITIES="OFF" -DXRAY="ON" ..
cmake --build . --config Release --target install

if %errorlevel%==0(
    echo penRed compilation completed
    cd ..\bindings\python\pyPenred
    pip install .
) else (
    echo "penRed compilation failed!"
    exit /b 1
)

if "%1"=="no-pause" (
    echo Install finished
) else (
    echo Install finished
    PAUSE
)
