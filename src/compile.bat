mkdir build

cd build
cmake -DWITH_DICOM="OFF" -DWITH_MULTI_THREADING="ON" -DWITH_MPI="OFF" -DWITH_LB="OFF" ..
cmake --build . --config Release --target install
cd ..

PAUSE
