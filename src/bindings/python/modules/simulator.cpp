//
//
//    Copyright (C) 2024 Universitat de València - UV
//    Copyright (C) 2024 Universitat Politècnica de València - UPV
//    Copyright (C) 2025-2026 Vicent Giménez Alventosa
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

#include "functions.hh"
#include "pen_simulation.hh"

#ifdef _PEN_XRAY_LIBS_
#include "x-ray.hh"
#endif

namespace py = pybind11;

using CompositionPair = std::pair<unsigned, double>;
using MaterialComposition = std::vector<CompositionPair>;
using MaterialData = std::pair<double, MaterialComposition>;
using MaterialList = std::vector<MaterialData>;
using FilterData = std::tuple<double, double, MaterialComposition>;
using FilterList = std::vector<FilterData>;

// + Tally results extraction

template<class TallyType, size_t I>
typename std::enable_if<I >= std::tuple_size<typename TallyType::ResultsTypes>::value, void>::type
tallyResults2numpy(const typename TallyType::ResultsTypes&, py::tuple&, const bool, const bool){}

template<class TallyType, size_t I>
typename std::enable_if<I < std::tuple_size<typename TallyType::ResultsTypes>::value, void>::type
tallyResults2numpy(const typename TallyType::ResultsTypes& r, py::tuple& pyRes, const bool extractInfo, const bool onlyEffective){
  pyRes[I] = result2numpy(std::get<I>(r), extractInfo, onlyEffective);
  tallyResults2numpy<TallyType, I+1>(r, pyRes, extractInfo, onlyEffective);
}

template<class TallyType>
py::tuple tallyExtractResults(const typename TallyType::ResultsTypes& r, const bool extractInfo, const bool onlyEffective){

  py::tuple pyRes(std::tuple_size<typename TallyType::ResultsTypes>::value);

  tallyResults2numpy<TallyType, 0>(r, pyRes, extractInfo, onlyEffective);

  return pyRes;
}

// + Tally search functions

template<size_t I = 0>
typename std::enable_if<I >= std::tuple_size<penred::tally::typesGenericTallies>::value, py::tuple>::type
getResults(const penred::tally::Results&, const std::string&, const bool, const bool){
  return py::tuple();
}  

template<size_t I = 0>
typename std::enable_if<I < std::tuple_size<penred::tally::typesGenericTallies>::value, py::tuple>::type
getResults(const penred::tally::Results& results, const std::string& tallyName, const bool extractInfo, const bool onlyEffective){

  using TallyType = typename std::tuple_element<I, penred::tally::typesGenericTallies>::type;

  const auto& tallyMap = results.read<I>();
  for(const auto& element : tallyMap){
    if(element.first.compare(tallyName) == 0){
      const typename TallyType::ResultsTypes& tallyResults = element.second;
      return tallyExtractResults<TallyType>(tallyResults, extractInfo, onlyEffective);
    }
  }
  return getResults<I+1>(results, tallyName, extractInfo, onlyEffective);
}

PYBIND11_MODULE(simulation,m){

  m.doc() = "penred simulation module";

  m.def("setConfigurationLog",
	[](const std::string& filename) -> void{
	  if(filename.empty())
	    penred::logs::logger::closeLogFile(penred::logs::CONFIGURATION);
	  else
	    penred::logs::logger::closeAndSetLogFile(penred::logs::CONFIGURATION, filename.c_str());
	},
	py::arg("filename"),
	R"(

Specify the file where the configuraiton logs will be redirected.

Args:
    filename (str) : Configuration log file name.
    
Returns:
    None

    )");
  
  m.def("setSimulationLog",
	[](const std::string& filename) -> void{
	  if(filename.empty())
	    penred::logs::logger::closeLogFile(penred::logs::SIMULATION);
	  else
	    penred::logs::logger::closeAndSetLogFile(penred::logs::SIMULATION, filename.c_str());
	},
	py::arg("filename"),
	R"(
Specify the file where the simulation logs will be redirected.

Args:
    filename (str) : Simulation log file name.
    
Returns:
    None

    )");

  m.def("errorMessage",
	[](const int& ierror) -> py::str{
	  return py::str(penred::simulation::errors::errorMessage(ierror));
	},
	py::arg("code"),
	R"(
Produces the error message associated to the specified error code.

Args:
    code (int) : Error code.
    
Returns:
    String containing the error message

    )");

  m.def("dict2SectionString", &dict2SectionString,
	py::arg("conf"),	
	R"(

Converts a dictionary to a compatible penRed configuration section string.

Args:
    conf (dict) : configuration dictionary to be converted.
    
Returns:
    String containing the converted dictionary

    )");

  m.def("configFile2YAML",
	[](const std::string& filename) -> py::str{
	  if(filename.empty())
	    return std::string("");
	  else{

	    //Parse configuration file
	    pen_parserSection config;
	    std::string errorLine;
	    unsigned long errorLineNum;
	    int err = parseFile(filename.c_str(),config,errorLine,errorLineNum);
	    
	    if(err != INTDATA_SUCCESS){
	      printf("Error parsing configuration.\n");
	      printf("Error code: %d\n",err);
	      printf("Error message: %s\n",pen_parserError(err));
	      printf("Error located at line %lu, at text: %s\n",
		     errorLineNum,errorLine.c_str());
	      return std::string("");
	    }

	    //Create YAML string
	    return config.stringifyYAML();
	  }
	},
	R"(

Reads a configuration file and returns a YAML string with the read information.

Args:
    filename (str) : File to be read.
    
Returns:
    String containing the information in YAML format

    )");

  py::class_<penred::simulation::simulator<pen_context>>(m, "simulator")
    .def(py::init<>())
    
    .def("configFromFile",
	 [](penred::simulation::simulator<pen_context>& obj,
	    const std::string& filename) -> void{
	   int err = obj.configFromFile(filename);
	   if(err != penred::simulation::errors::SUCCESS){
	     throw py::value_error(penred::simulation::errors::errorMessage(err));
	   }	   
	 },
	 py::arg("filename"),
	 R"(

Sets the entire configuration from a file matching the penRed internal data format.

Args:
    filename (str) : Configuration file.
    
Returns:
    None

Raises:
    ValueError: If file parsing fails.

    )")
    
    .def("configFromString",
	 [](penred::simulation::simulator<pen_context>& obj,
	    const std::string& configString) -> void{	   
	   int err = obj.configure(configString);
	   if(err != penred::simulation::errors::SUCCESS){
	     throw py::value_error(penred::simulation::errors::errorMessage(err));
	   }
	 },
	 py::arg("config"),
	 R"(

Sets the entire configuration from a string matching the penRed internal data format.

Args:
    config (str) : Text containing the configuration.
    
Returns:
    None

Raises:
    ValueError: If file parsing fails.

    )")
    
    .def("configure",
	 [](penred::simulation::simulator<pen_context>& obj,
	    const py::dict& config) -> void{
	   
	   std::string configString = dict2SectionString(config);	   
	   int err = obj.configure(configString);
	   if(err != penred::simulation::errors::SUCCESS){
	     throw py::value_error(penred::simulation::errors::errorMessage(err));
	   }	   
	 },
	 py::arg("config"),
	 R"(

Saves the entire simulation configuration from the provided dictionary.

Args:
    config (dict) : Dictionary with the simulation configuration.
    
Returns:
    None

Raises:
    ValueError: If configuraiton set fails.

    )")
    
    .def("configContext",
	 [](penred::simulation::simulator<pen_context>& obj,
	    const py::dict& config) -> void{

	   std::string configString = dict2SectionString(config);
	   
	   int err = obj.setContextConfig(configString);
	   if(err != penred::simulation::errors::SUCCESS){
	     throw py::value_error(penred::simulation::errors::errorMessage(err));
	   }	   
	 },
	 py::arg("config"),
	 R"(
Saves the context configuration from the provided dictionary.

Args:
    config (dict) : Dictionary with the context configuration.
    
Returns:
    None

Raises:
    ValueError: If configuraiton set fails.

    )")
    
    .def("configSource",
	 [](penred::simulation::simulator<pen_context>& obj,
	    const py::dict& config) -> void{

	   std::string configString = dict2SectionString(config);
	   
	   int err = obj.setSourceConfig(configString);
	   if(err != penred::simulation::errors::SUCCESS){
	     throw py::value_error(penred::simulation::errors::errorMessage(err));
	   }	   
	 },
	 py::arg("config"),
	 R"(

Saves the particle sources configuration from the provided dictionary.

Args:
    config (dict) : Dictionary with the source configuration.
    
Returns:
    None

Raises:
    ValueError: If configuraiton set fails.

    )")
    
    .def("configGeometry",
	 [](penred::simulation::simulator<pen_context>& obj,
	    const py::dict& config) -> void{

	   std::string configString = dict2SectionString(config);
	   
	   int err = obj.setGeometryConfig(configString);
	   if(err != penred::simulation::errors::SUCCESS){
	     throw py::value_error(penred::simulation::errors::errorMessage(err));
	   }
	 },
	 py::arg("config"),
	 R"(
Saves the geometry configuration from the provided dictionary.

Args:
    config (dict) : Dictionary with the geometry configuration.

Returns:
    None

Raises:
    ValueError: If configuraiton set fails.

    )")
    
    .def("configTally",
	 [](penred::simulation::simulator<pen_context>& obj,
	    const py::dict& config) -> void{

	   std::string configString = dict2SectionString(config);
	   
	   int err = obj.setTallyConfig(configString);
	   if(err != penred::simulation::errors::SUCCESS){
	     throw py::value_error(penred::simulation::errors::errorMessage(err));
	   }	   
	 },
	 py::arg("config"),
	 R"(

Saves the tallies configuration from the provided dictionary.

Args:
    config (dict) : Dictionary with the tallies configuration.

Returns:
    None

Raises:
    ValueError: If configuraiton set fails.

    )")
    
    .def("configVR",
	 [](penred::simulation::simulator<pen_context>& obj,
	    const py::dict& config) -> void{

	   std::string configString = dict2SectionString(config);
	   
	   int err = obj.setVRConfig(configString);
	   if(err != penred::simulation::errors::SUCCESS){
	     throw py::value_error(penred::simulation::errors::errorMessage(err));
	   }
	 },
	 py::arg("config"),
	 R"(

Saves the variance reduction configuration from the provided dictionary.

Args:
    config (dict) : Dictionary with the variance reduction configuration.

Returns:
    None

Raises:
    ValueError: If configuraiton set fails.

    )")
    
    .def("setSimConfig",
	 [](penred::simulation::simulator<pen_context>& obj,
	    const py::dict& config) -> void{

	   std::string configString = dict2SectionString(config);
	   
	   int err = obj.setSimConfig(configString);
	   if(err != penred::simulation::errors::SUCCESS){
	     throw py::value_error(penred::simulation::errors::errorMessage(err));
	   }	   
	 },
	 py::arg("config"),
	 R"(

Saves the simulation parameters from the provided dictionary.

Args:
    config (dict) : Dictionary with the variance simulation configuration.

Returns:
    None

Raises:
    ValueError: If configuraiton set fails.

    )")

    .def("addDumps",
         [](penred::simulation::simulator<pen_context>& obj, const std::vector<std::string>& filenames) -> void{

           int err;
           {
             //Release GIL for (probably) long simulation
             py::gil_scoped_release gil_release;
             err = obj.addDumps(filenames);
           }
           if(err != penred::simulation::errors::SUCCESS){
             throw py::value_error(penred::simulation::errors::errorMessage(err));
           }
         },
         R"(

Add dump files to generate either a unified results dump or ASCII files, depending on the configuration. The simulation object must already be configured with the same settings used to run the dumped simulations.

Args:
    filenames (list of str): List of dump filenames to add

Returns:
    None

Raises:
    ValueError: If invalid configurations are detected.

Example:
    .. code-block:: python

        #!/usr/bin/env python3
        import pyPenred

        sim = pyPenred.simulation.create()
        sim.configFromFile("config.in")

        sim.addDumps(["th0test.dump", "th1test.dump", "th3test.dump"])
    )")
    
    .def("simulate",
	 [](penred::simulation::simulator<pen_context>& obj, const bool async, const bool interactive) -> void {
	   if(async){
	     int err = obj.enableInteractive(interactive);
	     if(err != penred::simulation::errors::SUCCESS){
	       throw py::value_error(penred::simulation::errors::errorMessage(err));
	     }
	     err = obj.simulateAsync();
	     if(err != penred::simulation::errors::SUCCESS){
	       throw py::value_error(penred::simulation::errors::errorMessage(err));
	     }
	   }
	   else{
	     int err = obj.enableInteractive(false);
	     if(err != penred::simulation::errors::SUCCESS){
	       throw py::value_error(penred::simulation::errors::errorMessage(err));
	     }
	     {
	       //Release GIL for (probably) long simulation
	       py::gil_scoped_release gil_release;
	       err = obj.simulate();
	     }
	     if(err != penred::simulation::errors::SUCCESS){
	       throw py::value_error(penred::simulation::errors::errorMessage(err));
	     }
	   }
	 },
	 py::arg("async") = false,
	 py::arg("interactive") = false,
	 R"doc(
Executes the configured simulation in either blocking or non-blocking mode.

Args:
    async (bool, optional): Simulation execution mode.
        - False (default): Blocking mode (waits for completion)
        - True: Non-blocking mode (returns immediately)
    interactive (bool, optional): Simulation executes in interactive mode. Only available on asynchronous executions
        - False (default): The simulation starts and finish according to the configuration parameters.
        - True: Interactive. The simulation is configured and waits further instructions.

Returns:
    None

Raises:
    ValueError: If invalid configurations are detected pre-simulation.

Example:
    .. code-block:: python

        #!/usr/bin/env python3
        import pyPenred
        import time

        #Define configuration file
        configFile = "config.in"
    
        #Create simulation object
        simu = pyPenred.simulation.create()

        #Try to configure the simulation from the configuration file
        simu.configFromFile(configFile)
        
        print("Configuration set\n")

        #Start the simulation asynchronously
        simu.simulate(True)

        #Simulation started, check status every 20 seconds
        print("Simulation started\n")
        while simu.isSimulating():
            time.sleep(20)
            status = simu.stringifyStatus()
            for e in status:
                print(e)

        print("Simulation finished")

)doc")
    
    .def("simSpeeds",
	 [](penred::simulation::simulator<pen_context>& obj) -> py::tuple{
	   return py::tuple(py::cast(obj.simSpeeds()));
	 },
	 R"(

If a simulation is running, returns the simulation speeds, in histories per second, for each used thread

Args:
    None

Returns:
    tuple: A tuple containing the simulation speed of each thead

    )")
    
    .def("simulated",
	 [](penred::simulation::simulator<pen_context>& obj) -> py::tuple{
	   return py::tuple(py::cast(obj.simulated()));
	 },
	 R"(

If a simulation is running, returns the simulated histories for each thread at the current source.

Args:
    None

Returns:
    tuple: A tuple with the the pairs (simulated, to simulate) histories for each thread at the current source

)")
    .def("stringifyStatus",
	 [](penred::simulation::simulator<pen_context>& obj) -> py::tuple{
	   return py::tuple(py::cast(obj.stringifyStatus()));
	 },
	 R"(

If a simulation is running, returns a tuple with the stringified simulation status of each thread

Args:
    None

Returns:
    tuple: A tuple with a string for each thread describing its status.

)")
    
    .def("isSimulating",
	 [](penred::simulation::simulator<pen_context>& obj) -> py::bool_{
	   return obj.isSimulating();
	 },
	 R"(

Checks if a simulation is running.

Args:
    None

Returns:
    Bool: True if a simulation is running, false otherwise.

)")

    .def("forceFinish",
	 [](penred::simulation::simulator<pen_context>& obj) -> void{
	   int err = obj.forceFinish();
	   if(err != penred::simulation::errors::SUCCESS){
	     throw py::value_error(penred::simulation::errors::errorMessage(err));
	   }	   
	 },
	 R"(

Forces the simulation finish. Only works for non interactive simulations.

Args:
    None

Returns:
    None

Raises:
    ValueError: Simulation running in interactive mode.

)")
    .def("instructionEnd",
	 [](penred::simulation::simulator<pen_context>& obj) -> unsigned long long{
	   return obj.instructionEnd();
	 },
	 R"(

Appends a finish simulation instruction to the instructions queue. Once executed, the simulation is finished and the tally postprocessing starts.

Args:
    None

Returns:
    The numerical ID of the enqueued instruction.

)")
    .def("instructionSimulate",
	 [](penred::simulation::simulator<pen_context>& obj,
	    const std::string& sourceName, const double hists,
	    const std::vector<double>& translation,
	    const std::vector<double>& rotation) -> unsigned long long{

	   vector3D<double> t;
	   if(translation.size() >= 3){
	     t = vector3D<double>(translation[0], translation[1], translation[2]);
	   }

	   penred::Quaternion<double> r;
	   if(rotation.size() >= 4){
	     r = penred::Quaternion<double>(rotation[0],rotation[1],rotation[2],rotation[3]);
	   }
	   
	   return obj.instructionSimulate(sourceName, hists, t, r);
	 },
	 py::arg("source_name"),
	 py::arg("hists"),
	 py::arg("translation") = std::vector<double>(),
	 py::arg("quaternion") = std::vector<double>(),
	 R"(

Appends a simulate instruction to the instructions queue. Once executed, the specified source is simulated. If both, a rotation and a translation are provided, the former is applied first and the latter seconth. If provided, the three components must be specified.

Args:
    source_name (str): Name of the source to be simulated
    hists (float): Number of histories to simulate
    translation (list[float,float,float]): Specify the translation, in cm, to be applied to the original source position before the simulation starts. If no specified, no translation will be applied
    quaternion (list[float,float,float,float]): Specify the rotation quaternion to be applied to the original source before the simulation starts. If not provided, no extra rotation is applied to the original source.

Returns:
    The numerical ID of the enqueued instruction.

)")
    .def("getResults",
	 [](penred::simulation::simulator<pen_context>& obj,
	    const std::string& tallyName, const bool extractInfo, const bool onlyEffective) -> py::tuple{

	   py::tuple toRet = obj.processResults<py::tuple>
	     ([tallyName, extractInfo, onlyEffective](const penred::tally::Results& results){
	       return getResults<0>(results, tallyName, extractInfo, onlyEffective);
	     });
	   return toRet;
	 },
	 py::arg("tally_name"),
	 py::arg("extract_info") = false,
	 py::arg("only_effective") = false,
	 R"(

Gets the simulation results from the specified tally.

Args:
    tally_name (str): Name of the tally to get the results from.
    extract_info (bool): If enabled, the limits and dimensions information will be returned along with results values.
    only_effective (bool): If enabled, dimensions with a single bin are ignored, reducing the result's dimension.
Returns:
    On success, a tuple of numpy vectors storing the tally's specific results is returned. If the requested tally does not support retrieving results in that format, returned vectors will be empty.

)")    
    .def("instructionClear",
	 [](penred::simulation::simulator<pen_context>& obj) -> unsigned long long{
	   return obj.instructionClear();
	 },
	 R"(

Clears all pending instructions in queue.

Args:
    None

Returns:
    The last queued instruction numerical ID.

)")
    .def("instructionWait",
	 [](penred::simulation::simulator<pen_context>& obj,
	    const unsigned long long instructionID,
	    const unsigned long timeout) -> bool{
	   //Release GIL for waiting
	   py::gil_scoped_release gil_release;	   
	   return obj.instructionWait(instructionID, timeout);
	 },
	 py::arg("id"),
	 py::arg("timeout"),
	 R"(

Waits until the specified instruction ID is processed or the timeout is reached

Args:
    id (unsigned): ID of the instruction to wait.
    timeout (unsigned): Timeout in seconds.

Returns:
    Returns True if the instruction is processed before the timeout or False otherwise.

)")
    .def("instructionEndAndWait",
	 [](penred::simulation::simulator<pen_context>& obj,
	    const unsigned long timeout) -> bool{
	   //Release GIL for waiting
	   py::gil_scoped_release gil_release;	   
	   return obj.instructionEndAndWait(timeout);
	 },
	 py::arg("timeout"),
	 R"(

Enqueue an end instruction and waits until its completion or until reaching the timeout. In anycase, the end instruction sill enqueued until it is processed or cleared using a different call. 

Args:
    timeout (unsigned): Timeout in seconds.

Returns:
    Returns True if the end instruction is processed before the timeout or False otherwise.

)")
    .def("printResults",
	 [](penred::simulation::simulator<pen_context>& obj,
	    const std::string& prefix,
	    const unsigned nSigma,
	    const bool coordinates,
	    const bool binNumber,
	    const bool onlyEffective) -> void{
	   return obj.printResults(prefix.c_str(), nSigma, coordinates, binNumber, onlyEffective);
	 },
	 py::arg("prefix") = "",
	 py::arg("sigma") = 2,
	 py::arg("coordinates") = true,
	 py::arg("bins") = true,
	 py::arg("effective") = false,
	 R"(

Prints the latest available results to different files.

Args:
    prefix (str): Prefix applied to the filenames.
    sigma (unsigned): Number of printed sigmas for the uncertainties.
    coordinates (bool): Enable/disable coordinates printing.
    bins (bool): Enable/disable bin numbering printing.
    effective (bool): If enabled dimensions with a single bin will be ignored for printing.

Returns:
    Returns True if the end instruction is processed before the timeout or False otherwise.

)")
    .def("__repr__",
	 [](const penred::simulation::simulator<pen_context>& /*obj*/){
	   return std::string("<penred.simulator>");
	 })
    .doc() = R"doc(
A simulator class for penRed simulations.

This class handles particle transport and interaction modeling
in the penRed framework. It provides methods to configure and run
simulations, track particles, and retrieve results.
    )doc";

#ifdef _PEN_EMBEDDED_DATA_BASE_
  
  auto materials = m.def_submodule("materials");
  materials.doc() = "Materials utilities module";
  
  materials.def("create",
		[](const std::string& name, const double density,
		   MaterialComposition& comp,
		   const std::string& filename) -> void{

		  std::vector<penred::massFraction> composition;
		  for(const CompositionPair& e : comp){
		    composition.emplace_back(e.first, e.second);
		  }
		  std::string errorString;
		  int errMat;
		  {
		    //Release GIL during material creation
		    py::gil_scoped_release gil_release;
		    errMat = penred::penMaterialCreator::createMat(name, density, composition, errorString, filename);
		  }
		  if(errMat != 0)		  
		    throw py::value_error(errorString);
		},
		py::arg("name"),
		py::arg("density"),
		py::arg("composition"),
		py::arg("filename") = "",
		R"(
creates a material file based on the provided weight fraction composition. The composition should be given as a list of tuples with atomic number and weight fraction.

Args:
    name (str) : Name assigned to the material.
    density (float) : Material density in g/cm^3.
    composition (list of tuple) : Material composition list. Each element consists 
                                  of a 2D tuple with the corresponding element atomic number (Z) and weight fraction.
                                  For example, [(Z1, weight1), (Z2, weight2), ...].
    filename (str, optional) : File name for the generated material. If empty, it defaults to '${name}.mat'.

Returns:
    None

)");

  materials.def("createListed",
		[](const unsigned matID, const std::string& filename) -> void{
		  std::string errorString;
		  int err;
		  {
		    //Release GIL during material creation
		    py::gil_scoped_release gil_release;
		    err = penred::penMaterialCreator::createMat(matID, filename, errorString);
		  }
		  if(err != 0)
		    throw py::value_error(errorString);
		},
		py::arg("matID"),  // Material identifier
		py::arg("filename") = "",  // File name for the generated material, default empty string
		R"(

Creates a material file from the PENELOPE predefined material list.

Args:
    matID (int) : Material numerical identifier. The complete list can be found in the 'Predefined Materials'
                  appendix within the penred documentation.
    filename (str) : File name for the generated material.

Returns:
    None

)");
  
  
  materials.def("range",
		[](const std::vector<double>& energies,
		   const double density,
		   MaterialComposition& comp,
		   const unsigned verbose)->py::tuple{

		  // Create material
		  std::vector<penred::massFraction> composition;
		  for(const CompositionPair& e : comp){
		    composition.emplace_back(e.first, e.second);
		  }
		  std::string errorString;
		  int errMat;
		  {
		    //Release GIL during material creation
		    py::gil_scoped_release gil_release;
		    errMat = penred::penMaterialCreator::createMat("_range_", density, composition, errorString);
		  }
		  if(errMat != 0)
		    throw py::value_error(errorString);

		  // Create context
		  std::shared_ptr<pen_context> pcontext = createContext<pen_context>();

		  // Create context configuration
		  pen_parserSection contextConfig;
		  contextConfig.set("materials/mat1/number", 1);
		  contextConfig.set("materials/mat1/eabs/electron", 50.0);
		  contextConfig.set("materials/mat1/eabs/positron", 50.0);
		  contextConfig.set("materials/mat1/eabs/gamma", 50.0);
		  contextConfig.set("materials/mat1/C1", 0.2);
		  contextConfig.set("materials/mat1/C2", 0.2);
		  contextConfig.set("materials/mat1/WCC", 1.0e3);
		  contextConfig.set("materials/mat1/WCR", 1.0e3);
		  contextConfig.set("materials/mat1/filename", "_range_.mat");

		  // Context configuration
		  pen_parserSection matInfoSection;
		  if(pcontext->configure(1.0E9,
					contextConfig,
					matInfoSection,
					verbose) != pen_context::SUCCESS){
		    throw py::value_error("Error during configuration. Check context reports (context.rep)");
		  }

		  //Create the resulting python tuple
		  py::tuple results(constants::nParTypes);
		    
		  for(unsigned i = 0; i < constants::nParTypes; ++i)
		    {
		      std::vector<double> ranges(energies.size());
		      {
			//Release GIL during range calculations
			py::gil_scoped_release gil_release;
			for(unsigned j = 0; j < energies.size(); ++j)
			  {
			    ranges[j] = pcontext->range(energies[j], pen_KPAR(i), 0);
			  }
		      }
		      results[i] = py::tuple(py::cast(ranges));		      
		    }
	     
		  return results;
		},
		py::arg("energies"),
		py::arg("density"),
		py::arg("composition"),
		py::arg("verbose") = 1,
  	   R"doc(
Calculates the particle ranges (in cm) for a given material composition at specified energies.

Args:
    energies (double): List of energies, in eV, at which to calculate particle ranges.
    density (float): Material density in g/cm^3.
    composition (list of tuple) : Material composition list. Each element consists 
                                  of a 2D tuple with the corresponding element atomic number (Z) and weight fraction.
                                  For example, [(Z1, weight1), (Z2, weight2), ...].
    verbose (int, optional): Verbosity level.

Returns:
    tuple: A tuple containing the n tuples (one for each particle type), where each inner tuple contains the ranges, in cm, for all input energies. The order is: ((electron ranges), (gamma ranges), (positron ranges))

Example:
    .. code-block:: python

        # Run the simulation and store the results
        results = pyPenred.simulation.materials.range(energies=[100e3, 200e3, 50e3], density=1.0, composition=((1,0.112), (8,0.888)), verbose=1)

)doc");
  
 
  materials.def("rangeFromFile",
		[](const std::vector<double>& energies,
		   const std::string& filename,
		   const unsigned verbose)->py::tuple{
		   
		  // Create context
		  std::shared_ptr<pen_context> pcontext = createContext<pen_context>();

		  // Create context configuration
		  pen_parserSection contextConfig;
		  contextConfig.set("materials/mat1/number", 1);
		  contextConfig.set("materials/mat1/eabs/electron", 50.0);
		  contextConfig.set("materials/mat1/eabs/positron", 50.0);
		  contextConfig.set("materials/mat1/eabs/gamma", 50.0);
		  contextConfig.set("materials/mat1/C1", 0.2);
		  contextConfig.set("materials/mat1/C2", 0.2);
		  contextConfig.set("materials/mat1/WCC", 1.0e3);
		  contextConfig.set("materials/mat1/WCR", 1.0e3);
		  contextConfig.set("materials/mat1/filename",filename);

		  // Context configuration
		  pen_parserSection matInfoSection;
		  if(pcontext->configure(1.0E9,
					contextConfig,
					matInfoSection,
					verbose) != pen_context::SUCCESS){
		    throw py::value_error("Error during configuration. Check context reports (context.rep)");
		  }
		  
		  //Create the resulting python tuple
		  py::tuple results(constants::nParTypes);
		    
		  for(unsigned i = 0; i < constants::nParTypes; ++i)
		    {
		      std::vector<double> ranges(energies.size());
		      {
			//Release GIL during range calculations
			py::gil_scoped_release gil_release;			
			for(unsigned j = 0; j < energies.size(); ++j)
			  {
			    ranges[j] = pcontext->range(energies[j], pen_KPAR(i), 0);
			  }
		      }
		      results[i] = py::tuple(py::cast(ranges));		      
		    }
	     
		  return results;
		},
		py::arg("energies"),
		py::arg("filename"),
		py::arg("verbose") = 1,
  	   R"doc(
Calculates the particle ranges (in cm) for a given material file at specified energies.

Args:
    energies (double): List of energies, in eV, at which to calculate particle ranges.
    filename (str): Path to the material definition file.
    verbose (int, optional): Verbosity level.

Returns:
    tuple: A tuple containing the n tuples (one for each particle type), where each inner tuple contains the ranges, in cm, for all input energies. The order is: ((electron ranges), (gamma ranges), (positron ranges))

Example:
    .. code-block:: python

        # Run the simulation and store the results
        results = pyPenred.simulation.materials.rangeFromFile(energies=[100e3, 200e3, 50e3], filename=water.mat, verbose=1)

)doc");
  
  
 materials.def("mutrenInterval",
		[](const double density,
		   MaterialComposition& comp,
		   double emin,
		   double emax,
		   unsigned ebins, 
		   const double tolerance,
		   const double simTime)->py::tuple{
			   
		  // Create material
		  std::vector<penred::massFraction> composition;
		  for(const CompositionPair& e : comp){
		    composition.emplace_back(e.first, e.second);
		  }
		  std::string errorString;
		  if(penred::penMaterialCreator::createMat("_mutren_", density, composition, errorString) != 0)
		    throw py::value_error(errorString);
		 
		  if(emin<50.0)
			 emin=50.0;
			   
		  if(emax>1.0E9)
			 emax=1.0E9;
		
		  std::vector<double> muenData(ebins);
		  std::vector<double> EData(ebins);

		  {
		    //Release GIL during coefficient calculation
		    py::gil_scoped_release gil_release;		    
		    pen_muen::calculate(emin, emax, ebins, tolerance, simTime, "_mutren_.mat", EData, muenData);
		  }
		  
		  //Create the resulting python tuple
		  py::tuple results(2);
		  
		  results[0] = py::tuple(py::cast(EData));
		  results[1] = py::tuple(py::cast(muenData));
		 
		  return results;
		  		},
		py::arg("density"),
		py::arg("composition"),
		py::arg("emin") = 50.0,
		py::arg("emax") = 1.0E9,
		py::arg("ebins") = 100,
		py::arg("tolerance") = 0.1,
		py::arg("simTime") = 30,
  	   R"doc(
Calculates the mu_en coefficients for the specified material composition across a defined energy interval.

Args:
    density (float): Material density in g/cm^3.
    composition (list of tuple) : Material composition list. Each element consists 
                                  of a 2D tuple with the corresponding element atomic number (Z) and weight fraction.
                                  For example, [(Z1, weight1), (Z2, weight2), ...].
    emin (float, optional): Lower bound of energy range in eV. 
    emax (float, optional): Upper bound of energy range in eV. 
    ebins (unsigned, optional): Number of linear-spaced energy bins between emin and emax.
    tolerance (float, optional): Relative error to stop the simulation.
    simTime (float, optional): Allowed time (in seconds) to simulate each provided energy.
    
Returns:
    tuple: A tuple containing two tuples: The first with the energy bins and the second with the muen coefficients for the specified material.

Example:
    .. code-block:: python

        # Run the simulation and store the results
        results = pyPenred.simulation.materials.mutrenInterval(density=1.0, composition=((1,0.112), (8,0.888)), emin=100e3, emax=200e3, ebins=120, tolerance=0.1, simTime=30)

)doc");
			
 materials.def("mutren",
		[](const std::vector<double>& energies,
		   const double density,
		   MaterialComposition& comp,
		   const double tolerance,
		   const double simTime)->py::tuple{
			   
		  // Create material
		  std::vector<penred::massFraction> composition;
		  for(const CompositionPair& e : comp){
		    composition.emplace_back(e.first, e.second);
		  }
		  std::string errorString;
		  int errMat;
		  {
		    //Release GIL during material creation
		    py::gil_scoped_release gil_release;
		    errMat = penred::penMaterialCreator::createMat("_mutren_", density, composition, errorString);
		  }
		  if(errMat != 0)
		    throw py::value_error(errorString);
		
		  std::vector<double> muenData(energies.size());

		  {
		    //Release GIL during coefficient calculation
		    py::gil_scoped_release gil_release;
		    pen_muen::calculate(energies, tolerance, simTime, "_mutren_.mat", muenData);
		  }
		  
		  //Create the resulting python tuple
		  py::tuple results;
		  
		  results = py::tuple(py::cast(muenData));
		 
		  return results;
		},
		py::arg("energies"),
		py::arg("density"),
		py::arg("composition"),
		py::arg("tolerance") = 0.1,
		py::arg("simTime") = 30,
  	   R"doc(
Calculates the mu_en coefficients for the given material composition at specified energies.

Args:
    energies (list[float]): Array of energies, in eV.
    density (float): Material density in g/cm^3.
    composition (list of tuple) : Material composition list. Each element consists 
                                  of a 2D tuple with the corresponding element atomic number (Z) and weight fraction.
                                  For example, [(Z1, weight1), (Z2, weight2), ...].
    tolerance (float, optional): Relative error to stop the simulation.
    simTime (float, optional): Allowed time (in seconds) to simulate each provided energy.

Returns:
    tuple: A tuple with the muen coefficients for the specified material corresponding to the input energies.
	   For example, for a energy array provided with E1, E2, ..., En, returns:
	   (muen(E1), muen(E2), ..., muen(En))

Example:
    .. code-block:: python

        # Run the simulation and store the results
        results = pyPenred.simulation.materials.mutren(energies=muenE, density=1.0, composition=((1,0.112), (8,0.888)), tolerance=0.1, simTime=30)

)doc");
			
			
			  
 materials.def("mutrenFromFile",
		[](const std::vector<double>& energies, 
		   const std::string& filename,
		   const double tolerance,
		   const double simTime)->py::tuple{
		
		  std::vector<double> muenData(energies.size());
		
		  {
		    //Release GIL during coefficient calculation
		    py::gil_scoped_release gil_release;
		    pen_muen::calculate(energies, tolerance, simTime, filename.c_str(), muenData);
		  }
		  
		  //Create the resulting python tuple
		  py::tuple results;
		  
		  results = py::tuple(py::cast(muenData));
		 
		  return results;
		},
		py::arg("energies"),
		py::arg("filename"),
		py::arg("tolerance") = 0.1,
		py::arg("simTime") = 30,
  	   R"doc(
Calculates the mu_en coefficients for a given material file at specified energies.

Args:
    energies (list[float]): Array of energies, in eV.
    filename (str): Path to the material definition file.
    tolerance (float, optional): Relative error to stop the simulation.
    simTime (float, optional): Allowed time (in seconds) to simulate each provided energy.

Returns:
    tuple: A tuple with the muen coefficients for the specified material corresponding to the input energies.
           For example, for a energy array provided with E1, E2, ..., En, returns:
           (muen(E1), muen(E2), ..., muen(En))

Example:
    .. code-block:: python

        # Run the simulation and store the results
        results = pyPenred.simulation.materials.mutrenFromFile(energies=muenE, filename=water.mat, tolerance=0.1, simTime=30)

)doc");
  
#else
  // If embedded data base is not enabled, override submodule 'materials' to raise an exception on access
  m.def_submodule("materials").attr("__getattr__") =
    py::cpp_function([](py::object /* self */, const std::string& /* name */) -> py::object {
      throw py::import_error("Submodule 'materials' is disabled. "
			     "Recompile enabling EMBEDDED_DATA_BASE to include this module.");
      return py::none();
    });
#endif

#ifdef _PEN_XRAY_LIBS_
  
  auto xray = m.def_submodule("xray");
  xray.doc() = "x-ray utilities module";
    
  xray.def("anodeSim",
	   [](const double density, MaterialComposition& comp, const double anodeAngle,
	      const double kvp, const double nHists, const double minkeV,
	      const unsigned nBinsIn, const double pixelSize, const double maxAngle,
	      const bool saveDistributions, const unsigned verbose) -> py::tuple {

	     // Check parameters
	     if(nBinsIn == 0){
	       throw py::value_error("Error: Invalid number of bins: " + std::to_string(nBinsIn));
	     }

	     if(pixelSize <= 0.0){
	       throw py::value_error("Error: Invalid pixel size: " + std::to_string(pixelSize));
	     }

	     if(anodeAngle < 0.0 || anodeAngle >= 90.0){
	       throw py::value_error("Error: Valid anode angle values are [0,90)º. Invalid anode angle: " + std::to_string(anodeAngle));
	     }

	     const double beamEnergy = 1.0e3*kvp;
	     double minEnergy = 1.0e3*minkeV;
	     if(minEnergy < 0.0){
	       minEnergy = std::max(50.0, beamEnergy/10.0);
	     }
	     else if(minEnergy >= beamEnergy){
	       throw py::value_error("Error: Minimum energy mast be lesser than beam energy");
	     }

	     if(maxAngle <= 0.0 || maxAngle >= 90.0){
	       throw py::value_error("Error: Maximum angle must be within the range (0,90) degrees");
	     }

	     //Create anode material
	     std::vector<penred::massFraction> composition;
	     for(const CompositionPair& e : comp){
	       composition.emplace_back(e.first, e.second);
	     }
	     std::string errorString;

	     {
	       //Release GIL during material creation
	       py::gil_scoped_release gil_release;	       
	       penred::penMaterialCreator::createMat("anode", density, composition, errorString);
	     }

	     if(!errorString.empty()){
	       throw py::value_error("Error: Unable to create anode material: " + errorString);
	     }

	     //Run simulation
	     const unsigned long nBins = static_cast<unsigned long>(nBinsIn);
	     const unsigned long long nHist = static_cast<unsigned long long>(nHists);

	     //Create tallies
	     penred::measurements::measurement<double,1> spectrum;
	     penred::measurements::measurement<double,2> spatialDistrib;

	     spectrum.initFromLists({nBins}, {std::pair<double,double>(minEnergy,beamEnergy)});
  
	     //Simulate the anode
	     double dReg;
	     int err;
	     {
	       //Release GIL during anode simulation
	       py::gil_scoped_release gil_release;	       	       
	       err = penred::xray::simAnodeDistrib("anode.mat", beamEnergy, minEnergy, pixelSize, anodeAngle,
						   nHist, dReg, spectrum, spatialDistrib,
						   maxAngle, verbose);
	     }
	     if(err != penred::xray::errors::SUCCESS){
	       throw py::value_error(penred::xray::errors::message(err));	       
	     }

	     //Check if the results must be saved in disk
	     if(saveDistributions){
	       //Print the spectrum
	       FILE* fout = fopen("spectrum.dat", "w");
	       spectrum.print(fout, nHist, 2, true, false);
	       fclose(fout);

	       //Print the spatial distribution
	       fout = fopen("spatialDistrib.dat", "w");
	       spatialDistrib.print(fout, nHist, 2, true, true);
	       fclose(fout);	       
	     }
	     
	     //get results
	     penred::measurements::results<double,1> spectrumRes;
	     penred::measurements::results<double,2> spatialRes;
	     spectrum.results(nHist, spectrumRes);
	     spatialDistrib.results(nHist, spatialRes);
	     
	     py::array_t<double> npSpectrum =
	       py::array_t<double>(spectrumRes.readData().size(), spectrumRes.readData().data());
	     py::array_t<double> npSpectrumE =
	       py::array_t<double>(spectrumRes.readSigma().size(), spectrumRes.readSigma().data());
	     
	     py::array_t<double> npSpatialDistrib =
	       py::array_t<double>({
		   spatialRes.readDimBins()[1],
		   spatialRes.readDimBins()[0]
		 },
		 spatialRes.readData().data());
	     py::array_t<double> npSpatialDistribE =
	       py::array_t<double>({
		   spatialRes.readDimBins()[1],
		   spatialRes.readDimBins()[0]
		 },
		 spatialRes.readSigma().data());

	     //Get spectrum bin energy values
	     std::array<double, 2> eLimits =
	       {spectrum.readLimits()[0].first/1.0e3,
		spectrum.readLimits()[0].second/1.0e3};

	     py::array_t<double> npELimits =
	       py::array_t<double>(2, eLimits.data());
	     
	     //Get x,y bin values
	     std::vector<double> xLimits =
	       {spatialDistrib.readLimits()[0].first,
		spatialDistrib.readLimits()[0].second};
	     
	     std::vector<double> yLimits =
	       {spatialDistrib.readLimits()[1].first,
		spatialDistrib.readLimits()[1].second};

	     py::array_t<double> npXLimits =
	       py::array_t<double>(2, xLimits.data());
	     py::array_t<double> npYLimits =
	       py::array_t<double>(2, yLimits.data());
	     
	     
	     return py::make_tuple(npELimits, npSpectrum, npSpectrumE,
				   npXLimits, npYLimits,
				   npSpatialDistrib, npSpatialDistribE,
				   dReg);
	     
	   },
	   py::arg("density"),
	   py::arg("composition"),
	   py::arg("anode_angle"),
	   py::arg("kvp"),
	   py::arg("histories") = 1.0e4,
	   py::arg("min_energy") = -1.0,
	   py::arg("bins") = 200,
	   py::arg("pixel_size") = 0.1,
	   py::arg("max_angle") = 80.0,
	   py::arg("save_distributions") = true,
	   py::arg("verbose") = 1,
	   R"doc(
Simulates an electron beam impinging on an anode with no filter.

Args:
    density (float): Anode density in g/cm^3.
    composition (list[tuple[int, float]]): List of (Z, weight fraction) for material composition.
    anode_angle (float): Anode angle in degrees.
    kvp (float): Beam KVP value.
    histories (float, optional): Number of histories to simulate.
    min_energy (float, optional): Minimum energy to record in keV. Defaults to kvp/10.
    bins (int, optional): Number of spectrum bins.
    pixel_size (float, optional): Pixel size in cm.
    max_angle (float, optional): Maximum angle for scattered particles, in degrees.
    save_distributions(bool): If enabled, saves the spectrum and spatial distributions in files ready to be used by deviceSim.
    verbose (int, optional): Verbosity level.

Returns:
    tuple: A tuple with the following information, in order:
           - Energy limits for the generated spectrum (emin, emax) in keV
           - A 1D numpy array with the resulting gamma spectrum (prob/hist)
           - A 1D numpy array with the uncertainty for each gamma spectrum bin
           - X limits for spatial distribution (xmin, xmax) in cm
           - Y limits for spatial distribution (ymin, ymax) in cm
           - A 2D numpy array with the resulting gamma spatial distribution (prob/hist)
           - A 2D numpy array with the uncertainty for each spatial distribution bin
           - The distance between the anode collision point and the recorded spatial distribution center

Raises:
    ValueError: Incompatible value has been provided.

Example:
    .. code-block:: python

        #!/usr/bin/env python3
        import matplotlib.pyplot as plt
        import numpy as np
        import pyPenred

        # Run the simulation and store the results
        results = pyPenred.simulation.xray.anodeSim(density=10.0, composition=((74,1),), anode_angle=13, kvp=100.0, histories=1.0e5)

        # Extract energy and spatial ranges
        eLimits = results[0]
        xLimits = results[3]
        yLimits = results[4]

        # Produce X data for spectrum plot
        e_plot = np.linspace(eLimits[0], eLimits[1], len(results[1]))

        # Plot the spectrum
        plt.plot(e_plot, results[1])
        plt.title("Anode spectrum")
        plt.xlabel("keV")
        plt.ylabel("Prob/hist")
        plt.savefig('anode-spectrum.png')

        # Plot the spatial distribution
        plt.imshow(results[5],
                   extent=[xLimits[0], xLimits[1], yLimits[0], yLimits[1]],
                   aspect='auto',  
                   origin='lower', 
                   cmap='viridis')

        plt.colorbar()
        plt.xlabel("X cm")
        plt.ylabel("Y cm")
        plt.savefig('anode-spatial-distrib.png')

)doc");
  
  xray.def("deviceSim",
	   [](const std::array<double, 3>& detectorPos,
	      const double focalSpot, const double inherentFilterWidth,
	      const double minkeV, const double kvp,
	      const unsigned anodeZ, const double anodeAngle,
          const std::string& psfPath,
          const std::array<double, 3>& PSFTrans,
          const std::array<double, 3>& PSFRot,
          const std::string& spatialDistribPath,
          const std::string& energyDistribPath,
          const double distrib2detector,
	      const double source2filter,
	      const double source2detector,
	      const FilterList& filters,
	      const double nHists, const double maxTime,
	      const double dxDetector, double dyDetector, double dzDetector,
          const bool idealDetector, const MaterialData& detectorMat, 
	      const unsigned xbins, const unsigned ybins, const unsigned ebins,
	      const unsigned threads, const unsigned seedPair,
	      const py::dict& userTallies,
	      const py::dict& userGeometry,
	      const MaterialList& geoMats,
          const bool generateDetectorPSF,
          const bool generateFilteredPSF,
          const bool generateFilteredDistrib,
          const bool printConfig,
          const std::string& dump2read,
          const std::string& dump2write,
          const double dumpInterval,
	      const unsigned verbose) -> std::unique_ptr<penred::simulation::simulator<pen_context>> {
	     
	     // Check parameters
	     if(focalSpot < 0.0)
	       throw py::value_error("Error: Focal spot must be greater or equal to zero");

	     if(kvp < 1.0 || kvp > 1.0e3){
	       throw py::value_error("Error: Beam kvp must be within [1.0,1000]");
	     }
	     
	     double minEnergy = minkeV*1e3;
	     if(minEnergy < 0.0)
	       minEnergy = std::max(50.0, kvp/10.0);
	     else if(minEnergy < 50.0)
	       minEnergy = 50.0;
	     


	     if(minEnergy >= kvp*1.0e3)
	       throw py::value_error("Error: Beam kvp is lesser than minimum registered energy");
	     
	     if(anodeZ == 0)
	       throw py::value_error("Error: Invalid anode material");

	     if(anodeAngle <= 0.0 || anodeAngle >= 90.0)
	       throw py::value_error("Error: Anode angle value must be within (0,90) deg");

	     if(source2filter <= 0.0)
	       throw py::value_error("Error: Distance between source and first filter must be grater than zero");

	     if(source2detector <= 0.0)
	       throw py::value_error("Error: Distance between source and detector must be grater than zero");
         
	     if(!energyDistribPath.empty() && distrib2detector <= 0.0)
	       throw py::value_error("Error: Distance between distribution and detector must be grater than zero");
         
	     if(nHists <= 0.0)
	       throw py::value_error("Error: Number of histories must be grater than zero");

	     if(maxTime <= 0.0)
	       throw py::value_error("Error: Maximum simulation time must be greater than zero");

	     if(dxDetector <= 0.0 || dyDetector <= 0.0 || dzDetector <= 0.0)
	       throw py::value_error("Error: Detector size must be greater than zero");
         
	     if(xbins == 0 || ybins == 0 || ebins == 0){
	       throw py::value_error("Error: Number of bins must be greater than zero");
	     }

	     if(seedPair >= 1001)
	       throw py::value_error("Error: Seed pair index must be lesser or equal to 1000");

         if(generateDetectorPSF && generateFilteredPSF){
	       throw py::value_error("Error: Only one type of PSF can be generated per execution.");
         }
         
	     //Create config
	     pen_parserSection config;
	     config.set("simulation/detBins/nx", static_cast<int>(xbins));
	     config.set("simulation/detBins/ny", static_cast<int>(ybins));
	     config.set("simulation/eBins", static_cast<int>(ebins));
	     config.set("simulation/histories", nHists);
	     config.set("simulation/max-time", maxTime);
	     config.set("simulation/min-energy", minEnergy);
	     config.set("simulation/nthreads", static_cast<int>(threads));
	     config.set("simulation/seedPair", static_cast<int>(seedPair));

         //Dump config
         if(!dump2read.empty())
           config.set("simulation/dump/read", dump2read);
         if(!dump2write.empty())
           config.set("simulation/dump/write", dump2write);
         if(dumpInterval > 0.0)
           config.set("simulation/dump/time", dumpInterval);
         
	     if(psfPath.empty() && energyDistribPath.empty()){
	       config.set("simulation/sim-anode", true);
         }
         else if(!energyDistribPath.empty()){
	       config.set("simulation/sim-anode", false);
	       config.set("x-ray/source/distribution/spatial", spatialDistribPath);
	       config.set("x-ray/source/distribution/energy", energyDistribPath);
	       config.set("x-ray/source/distribution/distance", distrib2detector);
         }
	     else{
	       config.set("simulation/sim-anode", false);
	       config.set("x-ray/source/psf/path", psfPath);
           
           pen_parserArray PSFTransArray;
           PSFTransArray.append(PSFTrans[0]);
           PSFTransArray.append(PSFTrans[1]);
           PSFTransArray.append(PSFTrans[2]);
	       config.set("x-ray/source/psf/translation", PSFTransArray);

           pen_parserArray PSFRotArray;
           PSFRotArray.append(PSFRot[0]);
           PSFRotArray.append(PSFRot[1]);
           PSFRotArray.append(PSFRot[2]);
	       config.set("x-ray/source/psf/rotation", PSFRotArray);           
         }
	     
	     config.set("x-ray/anode/angle", anodeAngle);
	     config.set("x-ray/anode/z", static_cast<int>(anodeZ));
	     
	     config.set("x-ray/detector/dx", dxDetector);
	     config.set("x-ray/detector/dy", dyDetector);
	     config.set("x-ray/detector/dz", dzDetector);
	     config.set("x-ray/detector/ideal", idealDetector);
	     config.set("x-ray/detector/material", "detector.mat");
	     config.set("x-ray/distance/detector", source2detector);
	     config.set("x-ray/distance/filter", source2filter);
	     
	     config.set("x-ray/focal-spot", focalSpot);
	     config.set("x-ray/inherent-filter/width", inherentFilterWidth);
	     config.set("x-ray/kvp", kvp);

	     pen_parserArray detectorPosArray;
	     detectorPosArray.append(detectorPos[0]);
	     detectorPosArray.append(detectorPos[1]);
	     detectorPosArray.append(detectorPos[2]);
	     
	     config.set("x-ray/detector/position", detectorPosArray);

	     //Append filters
	     unsigned nextFilter = 1;
	     for(const FilterData& filter : filters){
	       std::string prefix = "x-ray/filters/filter" + std::to_string(nextFilter);
	       std::string matName = "filterMat" + std::to_string(nextFilter);
	       std::string matfilename = matName + ".mat";
	       config.set((prefix + "/mat-file").c_str(), matfilename);
	       config.set((prefix + "/width").c_str(), std::get<0>(filter));

	       //Create filter material
	       std::vector<penred::massFraction> composition;
	       for(const CompositionPair& e : std::get<2>(filter)){
		 composition.emplace_back(e.first, e.second);
	       }
	       std::string errorString;
	       int errMat;
	       {
		 //Release GIL during material creation
		 py::gil_scoped_release gil_release;
		 errMat = penred::penMaterialCreator::createMat(matName, std::get<1>(filter), composition,
								errorString, matfilename);
	       }
	       if(errMat != 0)		  
		 throw py::value_error(errorString);	       
	       ++nextFilter;
	     }

         //Create detector material, if needed
         if(idealDetector){
	       std::string errorString;
	       int errMat;
	       std::vector<penred::massFraction> composition;
	       for(const CompositionPair& e : detectorMat.second){
             composition.emplace_back(e.first, e.second);
	       }
           if(composition.empty()){
             //Release GIL during material creation
             py::gil_scoped_release gil_release;
             //Create Si by default
             errMat = penred::penMaterialCreator::createMat(14, "detector.mat", errorString);
           }
           else{
             //Release GIL during material creation
             py::gil_scoped_release gil_release;
             errMat = penred::penMaterialCreator::createMat("detector", detectorMat.first,
                                                            composition,
                                                            errorString,
                                                            "detector.mat");
	       }
	       if(errMat != 0)		  
             throw py::value_error(errorString);	       
         }

         //Append tally config
         if(!userTallies.empty()){
	       //Convert the dictionary to a configuration section
	       pen_parserSection tallyConfSec;
	       std::string errorString;
	       if(dict2section(userTallies, tallyConfSec, errorString) != 0)
             throw py::value_error(errorString);

	       if(config.addSubsection("tallies", tallyConfSec) != INTDATA_SUCCESS)
             throw py::value_error("Unable to set user tally configuration. "
                                   "Please, report this error");
         }

	     //Append geometry config
	     if(!userGeometry.empty()){
	       //Convert the dictionary to a configuration section
	       pen_parserSection geoConfSec;
	       std::string errorString;
	       if(dict2section(userGeometry, geoConfSec, errorString) != 0)
		 throw py::value_error(errorString);

	       if(config.addSubsection("geometry/config", geoConfSec) != INTDATA_SUCCESS)
		 throw py::value_error("Unable to set user geometry configuration. "
				       "Please, report this error");

	       //Append geometry materials
	       unsigned nextMat = 1;
	       for(const MaterialData& geoMat : geoMats){
		 std::string prefix = "geometry/materials/" + std::to_string(nextMat);
		 config.set((prefix + "/density").c_str(), geoMat.first);
		 config.set((prefix + "/number").c_str(), static_cast<int>(nextMat++));
		 prefix += "/elements/";
		 for(const CompositionPair& element : geoMat.second){
		   std::string key = prefix + std::to_string(element.first);
		   config.set(key, element.second);		 
		 }
	       }	       
	     }

         //Check PSF creation
         if(generateFilteredPSF){
           config.set("psf/filtered", true);
         }
         else if(generateDetectorPSF){
           config.set("psf/detected", true);
         }
         if(generateFilteredDistrib){
           config.set("distributions/filtered", true);
         }

	     if(printConfig){
	       FILE* fconf = fopen("simDevice.conf", "w");
	       fprintf(fconf, "%s\n", config.stringify().c_str());
	       fclose(fconf);
	     }

         //Set log files
         penred::logs::logger log;
         log.setConfigurationLogFile("config.log");
         log.setSimulationLogFile("simulation.log");

         //Create the simulation
         auto sim = std::make_unique<penred::simulation::simulator<pen_context>>();
         int err = penred::xray::constructSimDevice(config,*sim,verbose);  
	     if(err != penred::xray::errors::SUCCESS){
	       throw py::value_error(penred::xray::errors::message(err));	       
	     }
         
	     //Simulate the device
	     {
	       //Release GIL during device simulation
	       py::gil_scoped_release gil_release;
           err = sim->simulate();
	     }
	     if(err != penred::simulation::errors::SUCCESS){
	       throw py::value_error(penred::simulation::errors::errorMessage(err));
	     }

         return sim;	     
	   },
           py::arg("detector_position")       = std::array<double, 3>{0.0, 0.0, 0.0},
           py::arg("focal_spot")              = 0.0,
           py::arg("inherent_filter_width")   = -1.0,
           py::arg("min_energy")              = 10.0,
           py::arg("kvp")                     = 100.0,
           py::arg("anode_z")                 = 74,
           py::arg("anode_angle")             = 13.0,
           py::arg("psf_file")                = std::string(""),
           py::arg("psf_translation")         = std::array<double, 3>{0.0, 0.0, 0.0},
           py::arg("psf_rotation")            = std::array<double, 3>{0.0, 0.0, 0.0},
           py::arg("spatial_distribution")    = "",
           py::arg("energy_distribution")     = "",
           py::arg("distribution_to_detector")= -1.0,
           py::arg("source_to_filter")        = 7.0,
           py::arg("source_to_detector")      = 14.0,
           py::arg("filters")                 = FilterList{},
           py::arg("histories")               = 1.0e5,
           py::arg("max_time")                = 600.0,
           py::arg("detector_dx")             = 50.0,
           py::arg("detector_dy")             = 50.0,
           py::arg("detector_dz")             = 1.0,
           py::arg("detector_ideal")          = true,
           py::arg("detector_material")       = MaterialData{},
           py::arg("xbins")                   = 100,
           py::arg("ybins")                   = 100,
           py::arg("ebins")                   = 200,
           py::arg("threads")                 = 0,
           py::arg("seed_pair")               = 0,
           py::arg("user_tallies")            = py::dict(),
           py::arg("user_geometry")           = py::dict(),
           py::arg("geometry_materials")      = MaterialList{},
           py::arg("create_detected_psf")     = false,           
           py::arg("create_filtered_psf")     = false,           
           py::arg("create_filtered_distrib") = false,           
           py::arg("print_config")            = false,
           py::arg("dump_read")               = std::string(),
           py::arg("dump_write")              = std::string(),
           py::arg("dump_interval")           = -1.0,
           py::arg("verbose")                 = 1,
	   R"doc(
Simulates an electron beam impinging on an anode and records the resulting photon spectrum and spatial distribution. Alternatively, the simulation can be carried out from a Phase Space File (PSF) or sampling the initial photons from 2D (x,y) and 4D (E,x,y,z) spatial and energetic distributions respectively. Both, PSF and distributions can be obtained using this tool, enabling the 'create_filtered_psf' or 'create_filtered_distrib' options. 

Args:
    detector_position (tuple[float, float, float]): Detector position in cm. This point is interpreted as the center of the detector's top face in Z axis.
    focal_spot (float): Focal spot size of the beam in cm.
    inherent_filter_width (float): Width in cm for the inherent filter. Set to zero or a negative value to disable it.
    min_energy (float): Minimum photon energy to be recorded, in keV.
    kvp (float): Peak kilovoltage (kVp) of the beam.
    anode_z (int): Atomic number (Z) of the anode material.
    anode_angle (float): Angle of the anode in degrees.
    psf_file (str): Filename of the phase space file to be used as source.
    psf_translation (list): Translation to be applied to particles stored in the PSF (dx,dy,dz).
    psf_rotation (list): Euler angles to rotate the PSF's particles. Rotation is performed around the Z,Y,Z axis, in that order. Notice the detector position is below the source position (-Z direction).
    spatial_distribution (str): Path to the 2D spatial distribution (x,y) to be used as a source.
    energy_distribution (str): Path to the 4D energy distribution (e,x,y,z) to be used as a source.
    distribution_to_detector (float): Distance between the detector and the spatial source distribution
    source_to_filter (float): Distance from source to the first filter in cm.
    source_to_detector (float): Distance from source to the detector in cm.
    filters (list): List of additional filters beyond the inherent one. Each filter should be defined as:
        (width_cm, density_g_cm3, [(Z1, fraction1), (Z2, fraction2), ...])
    histories (float): Number of Monte Carlo histories (Beam electrons).
    max_time (float): Maximum simulation run time in seconds.
    detector_dx (float): Detector size in the X direction in cm.
    detector_dy (float): Detector size in the Y direction in cm.
    detector_dz (float): Detector depth in the Z direction in cm.
    detector_ideal (bool): Enables/disables ideal detection, i.e. perfect absorber.
    detector_material (list): Detector material for non-ideal detectors. The material must be defined as:
        (density, [(Z1, fraction1), (Z2, fraction2), ...])
    xbins (int): Number of bins (pixels) along the X axis.
    ybins (int): Number of bins (pixels) along the Y axis.
    ebins (int): Number of energy bins for the spectrum.
    threads (int): Number of threads to use. If 0, use all available CPU cores.
    seed_pair (int): Seed pair index for RNG initialization. Must be less than 1001.
    user_tallies (dict): Configuration of additional tallies provided by the user.
    user_geometry (dict): Configuration of additional geometry provided by the user.
    geometry_materials (list): List of extra materials used in the user-defined geometry. Each material should be defined as:
        (density, [(Z1, fraction1), (Z2, fraction2), ...])
    create_detected_psf (bool): If enabled, saves a PSF with the gammas reaching the detector. Only one PSF type can be created.
    create_filtered_psf (bool): If enabled, saves a PSF with the gammas just crossing the last filter. Only one PSF type can be created.
    create_filtered_distrib (bool): If enabled, saves the distribution of gammas just crossing the last filter. If no filter is configured (including inherent filtering), the scoring volume is set before the inherent filter position.
    print_config (bool): Prints the configuration using the penRed internal format. Only for debug purposes.
    dump_read (str): Specify the dumps' base-filename to be used to resume the simulation
    dump_write (str): Specify the dumps' base-filename to save the simulation status
    dump_interval (float): Specify the time interval, in seconds, to store simulation state dumps 
    verbose (int): Verbosity level (higher values provide more output).

Returns:
    Configured and executed 'simulator' object

Raises:
    ValueError: Incompatible value has been provided.

Example:
    .. code-block:: python

        import matplotlib.pyplot as plt
        import numpy as np
        import pyPenred

        sim = pyPenred.simulation.xray.deviceSim(ebins=100, inherent_filter_width=0.15, anode_angle=16, source_to_detector=100.0, max_time=600, histories=1.0e8)

        # Extract spatial distribution
        results = sim.getResults("SpatialDetector", extract_info=True, only_effective=True)

        # Extract spatial ranges
        xLimits = results[0][3]
        yLimits = results[0][2]

        # Plot it
        plt.figure(figsize=(6,5))
        plt.imshow(data,
                   origin='lower',
                   extent=[xLimits[0], xLimits[1], yLimits[0], yLimits[1]],
                   cmap='viridis',
                   aspect='equal')
        plt.xlabel('x (cm)')
        plt.ylabel('y (cm)')
        plt.colorbar(label='Total counts')
        plt.title('Energy-integrated detector image')
        plt.savefig("spatial-distrib.png", dpi=300)
        plt.close()

        # Extract mean detected spectrum
        resultsSpec = sim.getResults("SpectrumDetector", extract_info=True, only_effective=True)

        # Extract energy limits (in eV)
        eLimits = resultsSpec[0][2]
        e_plot = np.linspace(eLimits[0]/1000.0, eLimits[1]/1000.0, resultsSpec[0][0].shape[0])

        spectrum = resultsSpec[0][0]
        spectrumError = resultsSpec[0][1]

        plt.figure(figsize=(6,5))
        plt.errorbar(
            e_plot,
            spectrum,
            yerr=spectrumError,
            fmt='-',
            linewidth=1.2,
            ecolor='gray',
            elinewidth=0.8,
            capsize=2,
        )

        plt.xlabel("Energy (KeV)")
        plt.ylabel("Counts")
        plt.grid(True)
        plt.savefig("spectrum", dpi=300)
        plt.close()

)doc");
  
#else
  // If x-ray libriray is not enabled, override submodule 'xray' to raise an exception on access
  m.def_submodule("xray").attr("__getattr__") =
    py::cpp_function([](py::object /* self */, const std::string& /* name */) -> py::object {
      throw py::import_error("Submodule 'xray' is disabled. "
			     "Recompile enabling XRAY to include this module.");
      return py::none();
    });
#endif

}
