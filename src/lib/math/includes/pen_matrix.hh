
//
//
//    Copyright (C) 2019-2021 Universitat de València - UV
//    Copyright (C) 2019-2021 Universitat Politècnica de València - UPV
//
//    This file is part of PenRed: Parallel Engine for Radiation Energy Deposition.
//
//    PenRed is free software: you can redistribute it and/or modify
//    it under the terms of the GNU Affero General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    PenRed is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU Affero General Public License for more details.
//
//    You should have received a copy of the GNU Affero General Public License
//    along with PenRed.  If not, see <https://www.gnu.org/licenses/>. 
//
//    contact emails:
//
//        vicent.gimenez.alventosa@gmail.com
//        vicente.gimenez@uv.es
//    
//


#ifndef __PEN_MATRIX__
#define __PEN_MATRIX__

#include <stdexcept>
#include <cmath>

void matmul3D(const double A[9], double B[3]);
void matmul3D(const float A[9], float B[3]);

void createRotationZYZ(const double omega, const double theta, const double phi, double rotation[9]);

void rollAlign(const double u,
	       const double v,
	       const double w,
	       const double omega,
	       double rotation[9]);

void rollAlignf(const float u,
		const float v,
		const float w,
		const float omega,
		float rotation[9]);

#endif
