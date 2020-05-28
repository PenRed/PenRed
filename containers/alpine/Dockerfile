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

FROM alpine

WORKDIR /home

#Install, compile and remove all packages
RUN apk update && \
    apk add --no-cache git g++ make cmake && \
    git clone https://github.com/PenRed/PenRed.git && \
    mkdir PenRed/src/build && cd PenRed/src/build && \
    cmake -DWITH_DICOM="OFF" -DWITH_MULTI_THREADING="ON" -DWITH_NATIVE="OFF" ../ &&\
    make install && \
    cd /home && \
    cp PenRed/src/compiled/mains/penmain/pen_main /usr/bin/penred && \
    rm -rf PenRed && \
    apk del git make cmake g++ && \
    rm /var/cache/apk/*
#Install only the required libs
RUN apk add --no-cache libstdc++ libgcc
ENTRYPOINT ["/usr/bin/penred"]
