#
#
#    Copyright (C) 2020 Universitat de València - UV
#    Copyright (C) 2020 Universitat Politècnica de València - UPV
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
#        vicent.gimenez.alventosa@gmail.com (Vicent Giménez Alventosa)
#

# Download fedora as base image
FROM fedora:30

# Set root user
USER root

# Update image and install git
RUN dnf update -y
# Install git
RUN dnf install -y git
# Install compiler, cmake and make
RUN dnf install -y gcc-c++ make cmake
# Install DICOM reader library (dcmtk)
RUN dnf install -y dcmtk dcmtk-devel libjpeg libjpeg-devel libxml2 libxml2-devel openssl openssl-devel

#Clone PenRed repository and compile the code
WORKDIR root
RUN git clone https://github.com/PenRed/PenRed.git
RUN mkdir PenRed/src/build && cd PenRed/src/build && cmake -DWITH_DICOM="ON" -DWITH_MULTI_THREADING="ON" ../ && make install 

RUN cp PenRed/src/compiled/mains/penmain/pen_main /usr/bin/penred

# Create a PenRed user
RUN useradd -ms /bin/bash penred

# Change to penred user
USER penred
WORKDIR /home/penred

# Set entry point
ENTRYPOINT ["/usr/bin/penred"]
