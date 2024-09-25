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

std::string dict2SectionStringWithPrefix(const py::dict& dict, const std::string& prefixIn){

  //Convert a python dictionary to a string compatible with pen_parserSection
    
  std::string result;
  std::string prefix = prefixIn;
  if(!prefix.empty()){
    if(prefix.back() != '/'){
      //Append a slash
      prefix.append(1,'/');
    }
  }
    
  for(auto it : dict){
    //Check data type

    const std::string key = py::str(it.first).cast<std::string>();

    if(py::isinstance<py::bool_>(it.second)){
      //Boolean case
      if(it.second.cast<bool>())
	result += prefix + key + " true\n";
      else
	result += prefix + key + " false\n";
    }    
    else if(py::isinstance<py::int_>(it.second)){
      //Number case
      result += prefix + key + " " + std::to_string(it.second.cast<int>()) + "\n";
    }
    else if(py::isinstance<py::float_>(it.second)){
      //Float  case
      char aux[20];
      snprintf(aux, 20, " %15.5E\n", it.second.cast<double>());
      result += prefix + key + aux;
    }
    else if(py::isinstance<py::tuple>(it.second) ||
	    py::isinstance<py::list>(it.second)){
      //Array case
      result += prefix + key + " [";
      //Print each array element
      bool first = true;
      for(auto e : it.second){
	//Ensure the element is a number
	if(!py::isinstance<py::int_>(e) &&
	   !py::isinstance<py::float_>(e)){
	  printf("dict2SectionString: Error: The array at '%s' contains "
		 "a non numeric element. Will be skipped.\n",
		 (prefix + key).c_str());
	  continue;
	}
	if(py::isinstance<py::bool_>(e)){
	  if(it.second.cast<bool>())
	    result += (first ? " true" : ", true");
	  else
	    result += (first ? " false" : ", false");
	}
	else if(py::isinstance<py::int_>(e))
	  result += (first ? " " : ", ") + std::to_string(e.cast<int>());
	else{
	  char aux[25];
	  snprintf(aux, 25, "%s%15.5E\n", (first ? " " : ", "), e.cast<double>());	  
	  result += aux;
	}
	first = false;
      }
      result += " ]\n";
    }
    else if(py::isinstance<py::str>(it.second)){      
      //String case
      result += prefix + key + " \"" + it.second.cast<std::string>() + "\"\n";
    }
    else if(py::isinstance<py::dict>(it.second)){
      //Dictionary case. Add the prefix and parse it
      result += dict2SectionStringWithPrefix(it.second.cast<py::dict>(), prefix + key);
    }
    else{
      printf("dict2SectionString: Error: incompatible element at '%s'. "
	     "Only numbers, boolean, string, list, tuples and "
	     "dictionaries are allowed. Will be skipped.\n",
	     (prefix + key).c_str());
    }
  }
  return result;
}

std::string dict2SectionString(const py::dict& dict){
  return dict2SectionStringWithPrefix(dict, "");
}

PYBIND11_MODULE(simulation,m){

  m.doc() = "penred simulation module";

  m.def("setConfigurationLog",
	[](const std::string& filename){
	  if(filename.empty())
	    penred::logs::logger::setConfigurationLogFile(nullptr);
	  else
	    penred::logs::logger::setConfigurationLogFile(filename.c_str());
	}, "Sets configuration log file");
  m.def("setSimulationLog",
	[](const std::string& filename){
	  if(filename.empty())
	    penred::logs::logger::setSimulationLogFile(nullptr);
	  else
	    penred::logs::logger::setSimulationLogFile(filename.c_str());
	}, "Sets configuration log file");  

  m.def("dict2SectionString", &dict2SectionString, "Converts a dictionary to a compatible penRed configuration section string");

  py::class_<penred::simulation::simulator<pen_context>>(m, "simulator")
    .def(py::init<>())
    .def("configFromFile",
	 [](penred::simulation::simulator<pen_context>& obj,
	    const std::string& filename){
	   return obj.configFromFile(filename);
	 }, "Configure the whole simulation from a configuration file. The format must match the penred internal data structure.")
    .def("configFromString",
	 [](penred::simulation::simulator<pen_context>& obj,
	    const std::string& configString){	   
	   return obj.configure(configString);
	 }, "Configure the whole simulation from a string")
    .def("configure",
	 [](penred::simulation::simulator<pen_context>& obj,
	    const py::dict& config){
	   
	   std::string configString = dict2SectionString(config);	   
	   return obj.configure(configString);
	 }, "Configure the whole simulation using the provided dictionary")    
    .def("configContext",
	 [](penred::simulation::simulator<pen_context>& obj,
	    const py::dict& config){

	   std::string configString = dict2SectionString(config);
	   
	   return obj.setContextConfig(configString);
	 }, "Configure the context using the provided dictionary")
    .def("configSource",
	 [](penred::simulation::simulator<pen_context>& obj,
	    const py::dict& config){

	   std::string configString = dict2SectionString(config);
	   
	   return obj.setSourceConfig(configString);
	 }, "Configure particle sources using the provided dictionary")
    .def("configGeometry",
	 [](penred::simulation::simulator<pen_context>& obj,
	    const py::dict& config){

	   std::string configString = dict2SectionString(config);
	   
	   return obj.setGeometryConfig(configString);
	 }, "Configure the geometry using the provided dictionary")
    .def("configTally",
	 [](penred::simulation::simulator<pen_context>& obj,
	    const py::dict& config){

	   std::string configString = dict2SectionString(config);
	   
	   return obj.setTallyConfig(configString);
	 }, "Configure tallies using the provided dictionary")
    .def("configVR",
	 [](penred::simulation::simulator<pen_context>& obj,
	    const py::dict& config){

	   std::string configString = dict2SectionString(config);
	   
	   return obj.setVRConfig(configString);
	 })
    .def("setSimConfig",
	 [](penred::simulation::simulator<pen_context>& obj,
	    const py::dict& config){

	   std::string configString = dict2SectionString(config);
	   
	   return obj.setSimConfig(configString);
	 }, "Configure the simulation parameters using the provided dictionary")
    .def("simulate",
	 [](penred::simulation::simulator<pen_context>& obj, const bool async = false) -> int {
	   if(async){
	     //py::gil_scoped_release release;
	     return obj.simulateAsync();
	   }
	   else
	     return obj.simulate();
	 }, "Runs a simulation according to the provided configurations")
    .def("simSpeeds",
	 [](penred::simulation::simulator<pen_context>& obj){
	   return obj.simSpeeds();
	 }, "Returns a vector with the simulation speeds for each thread")
    .def("simulated",
	 [](penred::simulation::simulator<pen_context>& obj){
	   return obj.simulated();
	 }, "Returns a vector with the the pairs (simulated, to simulate) histories for each thread at the current source")
    .def("stringifyStatus",
	 [](penred::simulation::simulator<pen_context>& obj){
	   std::vector<std::string> r = obj.stringifyStatus();
	   auto resTuple = py::tuple(r.size());
	   for(size_t i = 0; i < r.size(); ++i)
	     resTuple[i] = py::str(r[i]);
	   return resTuple;
	 }, "Returns a vector with the stringified simulation state of each thread")        
    .def("isSimulating",
	 [](penred::simulation::simulator<pen_context>& obj){
	   return obj.isSimulating();
	 }, "Returns 'true' if the simulation is running and 'false' otherwise")
    .def("forceFinish",
	 [](penred::simulation::simulator<pen_context>& obj){
	   return obj.forceFinish();
	 }, "Forces the simulation finish")
    .def("__repr__",
	 [](const penred::simulation::simulator<pen_context>& /*obj*/){
	   return std::string("<penred.simulator>");
	 }); 
}
