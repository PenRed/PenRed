mkdir build
cd build
cmake -DWITH_DICOM="OFF" -DWITH_MULTI_THREADING="ON" -DWITH_MPI="OFF" -DWITH_LB="OFF" -DBUILD_PYTHON_MODULES="ON" -DBUILD_TESTS="OFF" -DBUILD_UTILITIES="OFF" ..
cmake --build . --config Release --target install

cd ..\bindings\python\pyPenred
pip install .

PAUSE
