
#
#
#    Copyright (C) 2024 Universitat de València - UV
#    Copyright (C) 2024 Universitat Politècnica de València - UPV
#    Copyright (C) 2024 Vicent Giménez Alventosa
#
#    This file is part of PenRed: Parallel Engine for Radiation Energy Deposition.
#
#    PenRed is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    PenRed is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with PenRed.  If not, see <https://www.gnu.org/licenses/>. 
#
#
#    contact emails:
#
#        vicent.gimenez.alventosa@gmail.com
#

rm -fr build
mkdir build
cd build
cmake -DWITH_DICOM="OFF" -DWITH_MPI="OFF" -DWITH_LB="OFF" -DBUILD_PYTHON_MODULES="ON" -DBUILD_TESTS="OFF" -DBUILD_UTILITIES="OFF" -DXRAY="ON" ..
cmake --build . --config Release --target install -j 4

if [ $? -eq 0 ]; then
    echo "penRed compilation completed"
    cd ../bindings/python/pyPenred
    pip install .
else
    echo "penRed compilation failed!"
fi
