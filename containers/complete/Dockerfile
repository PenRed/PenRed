#
#
#    Copyright (C) 2024 Universitat de València - UV
#    Copyright (C) 2024 Universitat Politècnica de València - UPV
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
FROM fedora:39

# Update image and install all dependencies, including DICOM support
RUN dnf update -y && dnf install -y git gcc-c++ make cmake dcmtk dcmtk-devel libjpeg libjpeg-devel libxml2 libxml2-devel openssl openssl-devel

#Clone PenRed repository, compile the code, save the main program and remove the build and source files
RUN cd && git clone https://github.com/PenRed/PenRed.git \
  && mkdir PenRed/src/build && cd PenRed/src/build \
  && cmake -DWITH_NATIVE="OFF" -DWITH_DICOM="ON" -DWITH_MULTI_THREADING="ON" ../ \
  && make install && cd \
  && cp PenRed/src/compiled/mains/penmain/pen_main /usr/bin/penred && rm -fr PenRed

# Add penred user and create simulation folder
RUN groupadd --gid 1000 penred \
  && useradd --uid 1000 --gid penred --create-home --no-log-init --shell /bin/bash penred \
  && mkdir /sim \
  && chown penred:penred /sim
USER penred

# Set working directory
WORKDIR /sim

# Set entry point
ENTRYPOINT ["/usr/bin/penred"]
