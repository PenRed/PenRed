
//
//
//    Copyright (C) 2019-2024 Universitat de València - UV
//    Copyright (C) 2019-2024 Universitat Politècnica de València - UPV
//    Copyright (C) 2025 Vicent Giménez Alventosa
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

 
#ifndef __PENRED_INTERNAL_DATA_CLASSES__
#define __PENRED_INTERNAL_DATA_CLASSES__

#include <utility>

#include "pen_parser.hh"
#include "pen_reader.hh"
#include "logger.hh"

//Constexpr string comparison
constexpr bool constexpr_strcmp(const char* a, const char* b) {
    return (*a == '\0' && *b == '\0') ? true :
           (*a == '\0' || *b == '\0') ? false :
           (*a != *b) ? false :
           constexpr_strcmp(a + 1, b + 1);
}

#endif
