//
//
//    Copyright (C) 2026 Vicent Giménez Alventosa
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
#include "geoHandler.hh"

namespace py = pybind11;


PYBIND11_MODULE(geometry,m){

  m.doc() = "penred geometry module";

  py::class_<penred::geometry::Handler>(m, "geometry")
    .def(py::init<>())

    .def("configFromFile",
	 [](penred::geometry::Handler& obj,
	    const std::string& filename,
	    const unsigned verbose,
	    const std::string& prefix) -> void{
	   
	   penred::errors::Error error = obj.configFromFile(filename.c_str(), verbose, prefix);
	   if(error){
	     throw py::value_error(error.stringify());
	   }
	   
	 },
	 py::arg("filename"),
	 py::arg("verbose") = 0,
	 py::arg("prefix") = "",
	 R"(

Configure the geometry from a file matching the penRed internal data format.

Args:
    filename (str) : Configuration file.
    verbose (int) : Verbose level.
    prefix (str) : Geometry section prefix within the configuration file.
    
Returns:
    None

Raises:
    ValueError: If configuration fails.

    )")
  
    .def("configFromString",
	 [](penred::geometry::Handler& obj,
	    const std::string& configString,
	    const unsigned verbose,
	    const std::string& prefix) -> void{
	   
	   penred::errors::Error error = obj.configure(configString, verbose, prefix);
	   if(error){
	     throw py::value_error(error.stringify());
	   }
	   
	 },
	 py::arg("config"),
	 py::arg("verbose") = 0,
	 py::arg("prefix") = "",
	 R"(

Configure the geometry from a string containing the configuration. This one must follows the penRed internal data format.

Args:
    config (str) : String containing the configuration. Must follows the penRed configuration format.
    verbose (int) : Verbose level.
    prefix (str) : Geometry section prefix within the configuration.
    
Returns:
    None

Raises:
    ValueError: If configuration fails.

    )")

    .def("configure",
	 [](penred::geometry::Handler& obj,
	    const py::dict& config,
	    const unsigned verbose) -> void{
	   
	   std::string configString = dict2SectionString(config);	   
	   penred::errors::Error error = obj.configure(configString, verbose);
	   if(error){
	     throw py::value_error(error.stringify());
	   }	   
	 },
	 py::arg("config"),
	 py::arg("verbose") = 0,
	 R"(

Configure the geometry from the provided dictionary

Args:
    config (dict) : Dictionary with the geometry configuration.
    verbose (int) : Verbose level.
    
Returns:
    None

Raises:
    ValueError: If configuraiton fails.

    )")

    .def("usedMat",
         [](penred::geometry::Handler& obj) -> py::list{

           //Read geometry
           std::shared_ptr<const wrapper_geometry> geo;
           penred::errors::Error error = obj.readGeometry(geo);
           if(error){
             throw py::value_error(error.stringify());
           }

           //Get used materials
           bool used[constants::MAXMAT+1];
           geo->usedMat(used);
           std::vector<unsigned> result;
           for(size_t i = 0; i < constants::MAXMAT+1; ++i){
             if(used[i]) result.push_back(i);
           }
           return py::cast(result);
         },
         R"(

Gets the materials used within the configured geometry

Args:
    None
    
Returns:
    A list with the index of each used material.

Raises:
    ValueError: If the geometry is not configured.

    )")    

    .def("numElements",
         [](penred::geometry::Handler& obj) -> unsigned{

           //Read geometry
           std::shared_ptr<const wrapper_geometry> geo;
           penred::errors::Error error = obj.readGeometry(geo);
           if(error){
             throw py::value_error(error.stringify());
           }

           return geo->getElements();
         },
         R"(

Return the total number of elements within the geoemtry. This value should match the
number of bodies in body-based geometries or the number of volume elements in mesh-based
geometries.

Args:
    None
    
Returns:
    The number of volume elements within the geometry.

Raises:
    ValueError: If the geometry is not configured.

    )")

    .def("elementsNDim",
         [](penred::geometry::Handler& obj) -> unsigned{

           //Read geometry
           std::shared_ptr<const wrapper_geometry> geo;
           penred::errors::Error error = obj.readGeometry(geo);
           if(error){
             throw py::value_error(error.stringify());
           }

           return geo->getElementsDim();
         },
         R"(

Return the number of dimensions used to organize the volumetric elements within the geoemtry. This value should be 1
 in body-based geometries. For example, for DICOM-based geometries this value is 3, because voxels are organizied in a regular
 3D mesh, corresponding to x, y and z axis.
geometries.

Args:
    None
    
Returns:
    The number of dimensions used to organize volume elements within the geometry.

Raises:
    ValueError: If the geometry is not configured.

    )")    

    .def("dimElements",
         [](penred::geometry::Handler& obj,
            const unsigned idim) -> unsigned{

           //Read geometry
           std::shared_ptr<const wrapper_geometry> geo;
           penred::errors::Error error = obj.readGeometry(geo);
           if(error){
             throw py::value_error(error.stringify());
           }

           return geo->getDimElements(idim);
         },
         R"(

Return the number of elements within the specified dimension.

Args:
    idim (int): Dimension index
    
Returns:
    The number of volume elements in the specified dimension.

Raises:
    ValueError: If the geometry is not configured.

    )")    

    .def("voxelSize",
         [](penred::geometry::Handler& obj) -> py::tuple{

           //Read geometry
           std::shared_ptr<const wrapper_geometry> geo;
           penred::errors::Error error = obj.readGeometry(geo);
           if(error){
             throw py::value_error(error.stringify());
           }

           //Check if is a voxelized-based geometry
           if(strcmp(geo->getType(), "VOXEL") == 0 || strcmp(geo->getType(), "DICOM") == 0){
             std::shared_ptr<const pen_voxelGeo> voxGeo = std::static_pointer_cast<const pen_voxelGeo>(geo);
             return py::make_tuple(voxGeo->xSize(), voxGeo->ySize(), voxGeo->zSize());
           }else{
             throw py::value_error("Configured geometry is not voxel-based");
           }
         },
         R"(

Return the voxel size in each dimension.

Args:
    None
    
Returns:
    A tuple with the voxel size, in cm, in each dimension [dx,dy,dz].

Raises:
    ValueError: If the geometry is not configured or if the geometry is not voxel-based.

    )")
    
    .def("numBodies",
         [](penred::geometry::Handler& obj) -> unsigned{

           //Read geometry
           std::shared_ptr<const wrapper_geometry> geo;
           penred::errors::Error error = obj.readGeometry(geo);
           if(error){
             throw py::value_error(error.stringify());
           }

           return geo->getBodies();
         },
         R"(

Return the geoemtry number of bodies

Args:
    None
    
Returns:
    The number of bodies within the geometry.

Raises:
    ValueError: If the geometry is not configured.

    )")    

    .def("bodyName",
         [](penred::geometry::Handler& obj,
            const unsigned ibody) -> std::string{

           //Read geometry
           std::shared_ptr<const wrapper_geometry> geo;
           penred::errors::Error error = obj.readGeometry(geo);
           if(error){
             throw py::value_error(error.stringify());
           }

           return geo->getBodyName(ibody);
         },
         py::arg("ibody"),
         R"(

Return the specified body's name

Args:
    ibody (int) : Body index within the geometry.
    
Returns:
    A string with the body name. The string is empty if the index is not used.

Raises:
    ValueError: If the geometry is not configured.

    )")

    .def("bodyIndex",
         [](penred::geometry::Handler& obj,
            const std::string& name) -> unsigned{

           //Read geometry
           std::shared_ptr<const wrapper_geometry> geo;
           penred::errors::Error error = obj.readGeometry(geo);
           if(error){
             throw py::value_error(error.stringify());
           }

           return geo->getIBody(name.c_str());
         },
         py::arg("name"),
         R"(

Return the specified body's index

Args:
    name (str) : Body name within the geometry.
    
Returns:
    The index of the specified body. The index is set to the total number of bodies if the name is not found.

Raises:
    ValueError: If the geometry is not configured.

    )")    
    
    .def("locate",
         [](penred::geometry::Handler& obj,
            const double x, const double y, const double z,
            const double time) -> py::tuple{

           //Read geometry
           std::shared_ptr<const wrapper_geometry> geo;
           penred::errors::Error error = obj.readGeometry(geo);
           if(error){
             throw py::value_error(error.stringify());
           }

           //Create the state
           pen_particleState state;
           state.PAGE = time;
           state.X = x;
           state.Y = y;
           state.Z = z;

           //Run the locate
           geo->locate(state);

           //Return material and body index
           return py::make_tuple(state.MAT, state.IBODY);
         },
         py::arg("x"),
         py::arg("y"),
         py::arg("z"),
         py::arg("time") = 0.0,
         R"(

Calculates the material and body index for the specified position and timestamp
within the configued geometry.

Args:
    x (float) : Position in the X axis, in cm.
    y (float) : Position in the Y axis, in cm.
    z (float) : Position in the Z axis, in cm.
    time (float) : Timestamp, in seconds (only used for animated geometries).
    
Returns:
    A tuple containing the material and body index, respectively, for the specified point.

Raises:
    ValueError: If the geometry is not configured.

    )")

    .def("step",
         [](penred::geometry::Handler& obj,
            const double x, const double y, const double z,
            const double u, const double v, const double w,
            const double time, const double ds) -> py::tuple{

           //Read geometry
           std::shared_ptr<const wrapper_geometry> geo;
           penred::errors::Error error = obj.readGeometry(geo);
           if(error){
             throw py::value_error(error.stringify());
           }

           //Normalize direction
           vector3D<double> dir(u,v,w);
           dir.normalize();

           //Create the state
           pen_particleState state;
           state.PAGE = time;
           state.X = x;
           state.Y = y;
           state.Z = z;
           state.U = dir.x;
           state.V = dir.y;
           state.W = dir.z;

           //First, locate the particle
           geo->locate(state);
           unsigned origMat = state.MAT;
           unsigned origBody = state.IBODY;

           //Next, move it
           double dsef, dstot;
           int ncross;
           geo->step(state, ds, dsef, dstot, ncross);

           //Return final material and body index, dsef, dstot and ncross
           return py::make_tuple(origMat, origBody, state.MAT, state.IBODY, dsef, dstot, ncross);
         },
         py::arg("x"),
         py::arg("y"),
         py::arg("z"),
         py::arg("u"),
         py::arg("v"),
         py::arg("w"),
         py::arg("time") = 0.0,
         py::arg("ds") = 1.0e35,
         R"(

Computes the specified movement through a previously configued geometry. The movement is halted if an
interaface is found.

Args:
    x,y,z (float) : Position within the geometry, in cm.
    u,v,w (float) : Movement's direction.
    time (float) : Timestamp, in seconds (only used for animated geometries).
    ds (float): Maximum distance to travel, in cm.
    
Returns:
    A tuple containing the original material and body index, the final material and body index, the distance traveled
    within the original material, the total traveled distance (including void regions), and the number of
    interfaces crossed, respectively.

Raises:
    ValueError: If the geometry is not configured.

    )")

    .def("renderX",
	 [](const penred::geometry::Handler& obj,
	    const float x, const float y, const float z,
        const float dy, const float dz,
        const unsigned ny, const unsigned nz,
        const float t,
        const unsigned nthreads) -> py::tuple{

       std::vector<unsigned int> renderMat;
       std::vector<unsigned int> renderBody;
	   penred::errors::Error error = obj.renderX(renderMat,renderBody,x,y,z,dy,dz,ny,nz,t,nthreads);
	   if(error){
	     throw py::value_error(error.stringify());
	   }

       py::array_t<unsigned> matArray = py::array_t<unsigned>({ny,nz}, renderMat.data());       
       py::array_t<unsigned> bodyArray = py::array_t<unsigned>({ny,nz}, renderBody.data());

       return py::make_tuple(matArray, bodyArray);
	 },
	 py::arg("x"),
	 py::arg("y"),
	 py::arg("z"),
	 py::arg("dy"),
	 py::arg("dz"),
	 py::arg("ny"),
	 py::arg("nz"),
	 py::arg("t") = 0.0,
	 py::arg("threads") = 1,
	 R"(

Renders the specified X plane from a previously configured geometry

Args:
    x (float): X coordinate of the render zone center, in cm
    y (float): Y coordinate of the render zone center, in cm
    z (float): Z coordinate of the render zone center, in cm
    dy (float): Pixel size in the Y axis, in cm
    dz (float): Pixel size in the Z axis, in cm
    ny (int): Number of pixels in the Y axis
    nz (int): Number of pixels in the Z axis
    t (float): Render time (only for animated geometries)
    threads (int): Number of threads used for rendering
    
Returns:
    A tuple containing two 2D numpy arrays, corresponding to the material and body assignment of each pixel, respectively.

Raises:
    ValueError: If geometry has not been configured previously.

    )")

    .def("renderY",
	 [](const penred::geometry::Handler& obj,
	    const float x, const float y, const float z,
        const float dx, const float dz,
        const unsigned nx, const unsigned nz,
        const float t,
        const unsigned nthreads) -> py::tuple{

       std::vector<unsigned int> renderMat;
       std::vector<unsigned int> renderBody;
	   penred::errors::Error error = obj.renderY(renderMat,renderBody,x,y,z,dx,dz,nx,nz,t,nthreads);
	   if(error){
	     throw py::value_error(error.stringify());
	   }

       py::array_t<unsigned> matArray = py::array_t<unsigned>({nx,nz}, renderMat.data());       
       py::array_t<unsigned> bodyArray = py::array_t<unsigned>({nx,nz}, renderBody.data());

       return py::make_tuple(matArray, bodyArray);
	 },
	 py::arg("x"),
	 py::arg("y"),
	 py::arg("z"),
	 py::arg("dx"),
	 py::arg("dz"),
	 py::arg("nx"),
	 py::arg("nz"),
	 py::arg("t") = 0.0,
	 py::arg("threads") = 1,
	 R"(

Renders the specified Y plane from a previously configured geometry

Args:
    x (float): X coordinate of the render zone center, in cm
    y (float): Y coordinate of the render zone center, in cm
    z (float): Z coordinate of the render zone center, in cm
    dx (float): Pixel size in the X axis, in cm
    dz (float): Pixel size in the Z axis, in cm
    nx (int): Number of pixels in the X axis
    nz (int): Number of pixels in the Z axis
    t (float): Render time (only for animated geometries)
    threads (int): Number of threads used for rendering
    
Returns:
    A tuple containing two 2D numpy arrays, corresponding to the material and body assignment of each pixel, respectively.

Raises:
    ValueError: If geometry has not been configured previously.

    )")

    .def("renderZ",
	 [](const penred::geometry::Handler& obj,
	    const float x, const float y, const float z,
        const float dx, const float dy,
        const unsigned nx, const unsigned ny,
        const float t,
        const unsigned nthreads) -> py::tuple{

       std::vector<unsigned int> renderMat;
       std::vector<unsigned int> renderBody;
	   penred::errors::Error error = obj.renderZ(renderMat,renderBody,x,y,z,dx,dy,nx,ny,t,nthreads);
	   if(error){
	     throw py::value_error(error.stringify());
	   }

       py::array_t<unsigned> matArray = py::array_t<unsigned>({nx,ny}, renderMat.data());       
       py::array_t<unsigned> bodyArray = py::array_t<unsigned>({nx,ny}, renderBody.data());

       return py::make_tuple(matArray, bodyArray);
	 },
	 py::arg("x"),
	 py::arg("y"),
	 py::arg("z"),
	 py::arg("dx"),
	 py::arg("dy"),
	 py::arg("nx"),
	 py::arg("ny"),
	 py::arg("t") = 0.0,
	 py::arg("threads") = 1,
	 R"(

Renders the specified Z plane from a previously configured geometry

Args:
    x (float): X coordinate of the render zone center, in cm
    y (float): Y coordinate of the render zone center, in cm
    z (float): Z coordinate of the render zone center, in cm
    dx (float): Pixel size in the X axis, in cm
    dz (float): Pixel size in the Z axis, in cm
    nx (int): Number of pixels in the X axis
    nz (int): Number of pixels in the Z axis
    t (float): Render time (only for animated geometries)
    threads (int): Number of threads used for rendering
    
Returns:
    A tuple containing two 2D numpy arrays, corresponding to the material and body assignment of each pixel, respectively.

Raises:
    ValueError: If geometry has not been configured previously.

    )")

    .def("toMesh",
	 [](const penred::geometry::Handler& obj,
	    const unsigned nx, const unsigned ny, const unsigned nz,
        const double dx, const double dy, const double dz,
        const double ox, const double oy, const double oz,
        const double t, const std::string& output,
        const unsigned granul,
        const unsigned smoothSteps,
        const double smoothFactor,
        const unsigned verbose) -> void{

	   penred::errors::Error error = obj.createMesh(nx, ny, nz, dx, dy, dz, ox, oy, oz, t,
                                                    output, granul, smoothSteps, smoothFactor, verbose);
	   if(error){
	     throw py::value_error(error.stringify());
	   }
	 },
         py::arg("nx"),
         py::arg("ny"),
         py::arg("nz"),
         py::arg("dx"),
         py::arg("dy"),
         py::arg("dz"),
         py::arg("ox"),
         py::arg("oy"),
         py::arg("oz"),
         py::arg("t") = 0.0,
         py::arg("output") = "geo",
         py::arg("granul") = 3,
         py::arg("smooth_steps") = 20,
         py::arg("smooth_factor") = 0.2,
         py::arg("verbose") = 0,
	 R"(

Creates a wavefront mesh scene from a previously configured geometry. Only the specified rectangular three-dimensional
area is used to generate the mesh. The area is defined by the origin and dimensions along each axis. The mesh is
assumed to extend in the positive direction of each axis. To construct the objects, the geometry is aproximated by an
auxiliary voxelized mesh, which parameters should be provided.

Args:
    nx (int): Number of auxiliary voxel planes used in the X axis
    ny (int): Number of auxiliary voxel planes used in the Y axis
    nz (int): Number of auxiliary voxel planes used in the Z axis
    dx (float): Auxiliary voxel size in the X axis, in cm
    dy (float): Auxiliary voxel size in the Y axis, in cm
    dz (float): Auxiliary voxel size in the Z axis, in cm
    ox (float): X coordinate of the mesh area origin, in cm
    oy (float): Y coordinate of the mesh area origin, in cm
    oz (float): Z coordinate of the mesh area origin, in cm
    t (float): Geometry time (only for animated geometries)
    output (str): Name used to create output files
    granul (int): Number of points, in each axis, used to determine the voxel object
    smooth_steps (int): Number of smooth steps performed to the final mesh
    smooth_factor (float): Strength factor for each smooth step
    
Returns:
    None. Creates two files, {output}.obj and {output}.mtl, containing the wavefront mesh and material data, respectively

Raises:
    ValueError: If geometry has not been configured previously.

    )");
    
}
