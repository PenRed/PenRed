//
//
//    Copyright (C) 2024 Universitat de València - UV
//    Copyright (C) 2024 Universitat Politècnica de València - UPV
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
//        vicent.gimenez.alventosa@gmail.com (Vicent Giménez Alventosa)
//        sanolgi@upvnet.upv.es (Sandra Oliver Gil)
//    
//

#include <pybind11/pybind11.h>

#include "pen_simulation.hh"

namespace py = pybind11;

PYBIND11_MODULE(simulator,m){

  m.doc() = "penred simulation module";

  py::class_<penred::simulation::simulator<pen_context>>(m, "simulator")
    .def(py::init<>())
    .def("configFromFile",
	 [](penred::simulation::simulator<pen_context>& obj,
	    const std::string& filename){
	   return obj.configFromFile(filename);
	 })
    .def("simulate",
	 [](penred::simulation::simulator<pen_context>& obj){
	   return obj.simulate();
	 })
    .def("__repr__",
	 [](const penred::simulation::simulator<pen_context>& /*obj*/){
	   return std::string("<penred.simulator>");
	 });
  
}
