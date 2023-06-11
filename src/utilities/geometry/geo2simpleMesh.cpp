
//
//
//    Copyright (C) 2023 Universitat de València - UV
//    Copyright (C) 2023 Universitat Politècnica de València - UPV
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
//    
//


#include <cstdio>
#include <vector>
#include <thread>
#include <atomic>
#include "palette.hh"
#include "pen_geometries.hh"
#include "MCtables.hh"

struct vertex{
  vector3D<float> v;
  std::vector<unsigned long> triangles;
  vector3D<float> normal;

  vertex(const vector3D<float>& vIn) : v(vIn){}
};

struct meshTriangle{
  vector3D<unsigned long> vi;
  vector3D<float> c;
  vector3D<float> normal;
};

int main(int argc, char** argv){

  if(argc < 2){
    printf("usage: %s config-filename\n",argv[0]);
    return 0;
  }

  unsigned verbose = 3;
  
  //**************************
  // Parse configuration
  //**************************

  //Parse configuration file
  pen_parserSection config;
  std::string errorLine;
  long unsigned errorLineNum;
  int err = parseFile(argv[1],config,errorLine,errorLineNum);
  
  if(err != INTDATA_SUCCESS){
    printf("Error parsing configuration.\n");
    printf("Error code: %d\n",err);
    printf("Error message: %s\n",pen_parserError(err));
    printf("Error located at line %lu, at text: %s\n",
	   errorLineNum,errorLine.c_str());
    return -1;
  }
  
  // Create origin geometry
  //************************

  //Get geometry section
  pen_parserSection geometrySection;
  err = config.readSubsection("geometry",geometrySection);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("Error: Configuration 'geometry' section doesn't exist.\n");
    }
    return -2;
  }

  //Get geometry type
  std::string geoType;
  if(geometrySection.read("type",geoType) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("createGeometry: Error: field 'geometry/type' not specified. String expected.\n");
    }
    return -3;
  }  
  
  wrapper_geometry* geometry = nullptr;
  //Create geometry
  geometry = penGeoRegister_create(geoType.c_str());
  if(geometry == nullptr){
    if(verbose > 0){
      printf("Error creating geometry instance of type '%s'\n", geoType.c_str());
    }
    return -4;
  }

  //Configure geometry  
  geometry->name.assign("geometry");    
  err = geometry->configure(geometrySection,verbose);
  if(err != 0){
    printf("Error: Unable to configure geometry\n");
    return -5;
  }

  //Get number of bodies
  unsigned long nBodies = geometry->getBodies();
  //For quadric geometries, remove the last one
  if(geoType.compare("PEN_QUADRIC") == 0){
    nBodies = nBodies > 0 ? nBodies-1 : 0;
  }

  if(nBodies == 0){
    printf("Error: Geometry does not contain any body\n");
    return -6;
  }
  
  //Create containers for body meshes
  std::vector<std::vector<vertex>> bodyVertex(nBodies);
  std::vector<std::vector<meshTriangle>> bodyTriangles(nBodies);
  
  //Check errors
  if(geometry->configureStatus() != 0){
    if(verbose > 0)
      printf("Error: Fail on geometry configuration.\n");
    return -5;
  }

  // Create voxelized geometry
  //***************************

  //Read voxelized geometry parameters
  int nx,ny,nz;
  double dx,dy,dz;
  double ox,oy,oz;
  int granul;
  const unsigned maxGranul = 200;
  
  //Number of voxels
  if(config.read("nx",nx) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("createGeometry: Error: field 'voxelized/nx' not specified. integer expected.\n");
    }
    return -6;
  }
  if(config.read("ny",ny) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("createGeometry: Error: field 'voxelized/ny' not specified. integer expected.\n");
    }
    return -6;
  }
  if(config.read("nz",nz) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("createGeometry: Error: field 'voxelized/nz' not specified. integer expected.\n");
    }
    return -6;
  }

  //Voxel sizes
  if(config.read("dx",dx) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("createGeometry: Error: field 'voxelized/dx' not specified. Double expected.\n");
    }
    return -7;
  }
  if(config.read("dy",dy) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("createGeometry: Error: field 'voxelized/dy' not specified. Double expected.\n");
    }
    return -7;
  }
  if(config.read("dz",dz) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("createGeometry: Error: field 'voxelized/dz' not specified. Double expected.\n");
    }
    return -7;
  }

  //Origin
  if(config.read("ox",ox) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("createGeometry: Error: field 'voxelized/ox' not specified. Double expected.\n");
    }
    return -8;
  }
  if(config.read("oy",oy) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("createGeometry: Error: field 'voxelized/oy' not specified. Double expected.\n");
    }
    return -8;
  }
  if(config.read("oz",oz) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("createGeometry: Error: field 'voxelized/oz' not specified. Double expected.\n");
    }
    return -8;
  }  
  

  //Read granularity
  if(config.read("granularity",granul) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("createGeometry: Error: field 'voxelized/granularity' not specified. integer expected.\n");
    }
    return -9;
  }


  //Read smooth parameters
  unsigned smoothSteps = 50;
  int smoothStepsi;  
  if(config.read("smooth/steps",smoothStepsi) == INTDATA_SUCCESS){
    if(smoothStepsi >= 0)
      smoothSteps = static_cast<unsigned>(smoothStepsi);
  }

  float smoothFactor = 0.2;
  double smoothFactord;  
  if(config.read("smooth/factor",smoothFactord) == INTDATA_SUCCESS){
    if(smoothFactord > 0)
      smoothFactor = static_cast<float>(smoothFactord);
  }
  
  //Check values
  if( nx < 1 || ny < 1 || nz < 1 ){
    printf("Error: Number of voxels must be, as least, 1 on each axis.\n");
    printf("            nx: %d\n",nx);
    printf("            ny: %d\n",ny);
    printf("            nz: %d\n",nz);
    return -9;
  }

  if( dx < 1.0e-15 || dy < 1.0e-15 || dz < 1.0e-15 ){
    printf("Error: Voxel sizes must be greater than zero.\n");
    printf("            dx: %12.4E\n",dx);
    printf("            dy: %12.4E\n",dy);
    printf("            dz: %12.4E\n",dz);
    return -10;
  }

  if(granul <= 0 || granul > (int)maxGranul){
    printf("Error: Granularity must be, as least, 1 and less than %u\n",maxGranul);
    printf("   granularity: %d\n",granul);
    return -11;
  }
  
  //Create voxelized geometry
  unsigned long* voxBodies = nullptr;
  
  //Allocate memory for two slices
  voxBodies = new unsigned long[nx*ny*2];

  if(voxBodies == nullptr){
    printf("Error: Unable to allocate memory for voxel data\n");
    return -12;
  }

  //Fill voxel matrix
  double subSteps[3] = {dx/(double)granul,dy/(double)granul,dz/(double)granul};

  // ** Set thread variables
  unsigned int nConcurrency = std::max(std::thread::hardware_concurrency(), 2u);
  std::vector<std::thread> threads;
  threads.reserve(nConcurrency);
  std::atomic<unsigned long> atomicCounter(0);
  //****************************
  
  std::vector<unsigned long> firstVertexPrevZPlane(nBodies,0);
  std::vector<unsigned long> firstVertexPrev2ZPlane(nBodies,0);
  
  //Fill the first plane with multiple threads
  for(size_t ith = 0; ith < nConcurrency; ++ith){

    threads.push_back(std::thread([=,&atomicCounter]{

      int j = atomicCounter++;
      while(j < ny){
	for(int i = 0; i < nx; ++i){

	  if(granul > 1){
	    unsigned pointsPerBody[pen_geoconst::NB];
	    for(unsigned ibody = 0; ibody < pen_geoconst::NB; ++ibody)
	      pointsPerBody[ibody] = 0;

	    //Set initial state position
	    pen_particleState state;	    
	    double xInit = (double)i*dx-subSteps[0]*0.5+ox;
	    double yInit = (double)j*dy-subSteps[1]*0.5+oy;
	    double zInit = -subSteps[2]*0.5+oz;
	    state.X = xInit;
	    state.Y = yInit;
	    state.Z = zInit;
	    for(int k2 = 0; k2 < granul; ++k2){
	      state.Z += subSteps[2];
	      state.Y = yInit;
	      for(int j2 = 0; j2 < granul; ++j2){
		state.Y += subSteps[1];
		state.X = xInit;
		for(int i2 = 0; i2 < granul; ++i2){
		  state.X += subSteps[0];

		  //Locate the state
		  geometry->locate(state);
		  ++pointsPerBody[state.IBODY];
		}
	      }
	    }

	    //Set material
	    unsigned maxPoints = 0;
	    unsigned long dominantBody = nBodies;
	    for(unsigned ibody = 0; ibody < nBodies; ++ibody){
	      if(pointsPerBody[ibody] > maxPoints){
		maxPoints = pointsPerBody[ibody];
		dominantBody = ibody;
	      }
	    }
	    voxBodies[j*nx + i] = dominantBody;
	
	  }else{

	    //Set auxiliary state position
	    pen_particleState state;
	    state.X = ((double)i)*dx+ox;
	    state.Y = ((double)j)*dy+oy;
	    state.Z = oz;
	  
	    //Locate the state
	    geometry->locate(state);
	    voxBodies[j*nx + i] = state.IBODY;
	  }
	}
	j = atomicCounter++;
      }
    }));
  }

  //Join threads
  for(size_t ith = 0; ith < nConcurrency; ++ith){
    threads[ith].join();
  }
  //Clear thread vector
  threads.clear();

  //Calculate the mesh
  unsigned long* firstPlane = voxBodies;  
  unsigned long* lastPlane = &voxBodies[nx*ny];

  printf("\n");

  for(int k = 1; k < nz; ++k){

    printf("\rProcessed planes:%04d/%04d",k,nz);
    
    //Init atomic counter
    atomicCounter = 0;
    
    //Fill the next plane with multiple threads
    for(size_t ith = 0; ith < nConcurrency; ++ith){
      
      threads.push_back(std::thread([=,&atomicCounter]{

	//Get next plane
	int j = atomicCounter++;
	while(j < ny){
	  
	  for(int i = 0; i < nx; ++i){

	    if(granul > 1){
	      unsigned pointsPerBody[pen_geoconst::NB];
	      for(unsigned ibody = 0; ibody < pen_geoconst::NB; ++ibody)
		pointsPerBody[ibody] = 0;

	      //Set initial state position
	      pen_particleState state;	    
	      double xInit = (double)i*dx-subSteps[0]*0.5+ox;
	      double yInit = (double)j*dy-subSteps[1]*0.5+oy;
	      double zInit = (double)k*dz-subSteps[2]*0.5+oz;
	      state.X = xInit;
	      state.Y = yInit;
	      state.Z = zInit;
	      for(int k2 = 0; k2 < granul; ++k2){
		state.Z += subSteps[2];
		state.Y = yInit;
		for(int j2 = 0; j2 < granul; ++j2){
		  state.Y += subSteps[1];
		  state.X = xInit;
		  for(int i2 = 0; i2 < granul; ++i2){
		    state.X += subSteps[0];

		    //Locate the state
		    geometry->locate(state);
		    ++pointsPerBody[state.IBODY];
		  }
		}
	      }

	      //Set material
	      unsigned maxPoints = 0;
	      unsigned long dominantBody = nBodies;
	      for(unsigned ibody = 0; ibody < nBodies; ++ibody){
		if(pointsPerBody[ibody] > maxPoints){
		  maxPoints = pointsPerBody[ibody];
		  dominantBody = ibody;
		}
	      }
	      lastPlane[j*nx + i] = dominantBody;	  
	    }else{

	      //Set auxiliary state position
	      pen_particleState state;	    
	      state.X = ((double)i)*dx+ox;
	      state.Y = ((double)j)*dy+oy;
	      state.Z = ((double)k)*dz+oz;
	  
	      //Locate the state
	      geometry->locate(state);
	      lastPlane[j*nx + i] = state.IBODY;
	    }
	  }
	  j = atomicCounter++;
	}

      }));
    }

    //Join threads
    for(size_t ith = 0; ith < nConcurrency; ++ith){
      threads[ith].join();
    }
    //Clear thread vector
    threads.clear();
    

    //Store initial number of vertex until this plane 
    for(unsigned ibody = 0; ibody < nBodies; ++ibody){
      firstVertexPrevZPlane[ibody] = bodyVertex[ibody].size();
    }
    

    float z1  = static_cast<float>(k-1)*dz + oz;
    float z15  = (static_cast<float>(k)-0.5)*dz + oz;
    float z2  = static_cast<float>(k)*dz + oz;
      
    //Set vertex, avoiding borders
    for(int j = 0; j < ny-1; ++j){

      float y1  = static_cast<float>(j)*dy + oy;
      float y15  = (static_cast<float>(j)+0.5)*dy + oy;
      float y2  = static_cast<float>(j+1)*dy + oy;

      for(int i = 0; i < nx-1; ++i){

	float x1  = static_cast<float>(i)*dx + ox;
	float x15  = (static_cast<float>(i)+0.5)*dx + ox;
	float x2  = static_cast<float>(i+1)*dx + ox;

	//Get vertex references values
	const unsigned long v1 = firstPlane[j*nx + i];
	const unsigned long v2 = firstPlane[j*nx + i + 1];
	const unsigned long v3 = firstPlane[(j+1)*nx + i + 1];
	const unsigned long v4 = firstPlane[(j+1)*nx + i];
	const unsigned long v5 = lastPlane[j*nx + i];
	const unsigned long v6 = lastPlane[j*nx + i + 1];
	const unsigned long v7 = lastPlane[(j+1)*nx + i + 1];
	const unsigned long v8 = lastPlane[(j+1)*nx + i];
	
	//Construct possible triangle vertex
	vector3D<float> e[12];

	//Front
	e[0].x = x15; e[0].y = y1; e[0].z = z1;

	e[1].x = x2; e[1].y = y15; e[1].z = z1;

	e[2].x = x15; e[2].y = y2; e[2].z = z1;

	e[3].x = x1; e[3].y = y15; e[3].z = z1;

	//Back
	e[4].x = x15; e[4].y = y1; e[4].z = z2;

	e[5].x = x2; e[5].y = y15; e[5].z = z2;

	e[6].x = x15; e[6].y = y2; e[6].z = z2;

	e[7].x = x1; e[7].y = y15; e[7].z = z2;

	//Center
	e[8].x = x1; e[8].y = y1; e[8].z = z15;

	e[9].x = x2; e[9].y = y1; e[9].z = z15;

	e[10].x = x2; e[10].y = y2; e[10].z = z15;

	e[11].x = x1; e[11].y = y2; e[11].z = z15;
	  
	//Iterate over geometry bodies
	for(unsigned ibody = 0; ibody < nBodies; ++ibody){

	  //Get cube index
	  int cubeIndex = 0;
	  if(v1 == ibody) cubeIndex |= (1 << 0);
	  if(v2 == ibody) cubeIndex |= (1 << 1);
	  if(v3 == ibody) cubeIndex |= (1 << 2);
	  if(v4 == ibody) cubeIndex |= (1 << 3);
	  if(v5 == ibody) cubeIndex |= (1 << 4);
	  if(v6 == ibody) cubeIndex |= (1 << 5);
	  if(v7 == ibody) cubeIndex |= (1 << 6);
	  if(v8 == ibody) cubeIndex |= (1 << 7);

	  //Check if it is completely inside or outside
	  if(edgeTable[cubeIndex] == 0) continue;
	  
	  //Iterate over cube triangles according to its pattern
	  for(int n = 0; triTable[cubeIndex][n] != -1; n+=3){
	    unsigned long nVertex = bodyVertex[ibody].size();

	    const vector3D<float> center = (e[n] + e[n+1] + e[n+2])/3.0f;

	    //Check if some vertex already exists to avoid repetitions
	    vector3D<long int> triVertex = {-1,-1,-1};

	    //V0
	    for(size_t iv = firstVertexPrev2ZPlane[ibody]; iv < nVertex; ++iv){
	      if(bodyVertex[ibody][iv].v.dist(e[triTable[cubeIndex][n+2]]) < 1.0e-3){
		triVertex.x = iv;
		//Add the triangle to this vertex
		bodyVertex[ibody][iv].triangles.push_back(bodyTriangles[ibody].size());
		break;
	      }
	    }

	    if(triVertex.x < 0){
	      triVertex.x = nVertex++;
	      bodyVertex[ibody].push_back(e[triTable[cubeIndex][n+2]]);  //v0
	      //Add the triangle to this vertex
	      bodyVertex[ibody].back().triangles.push_back(bodyTriangles[ibody].size());
	    }

	    
	    //V1
	    for(size_t iv = firstVertexPrev2ZPlane[ibody]; iv < nVertex; ++iv){
	      if(bodyVertex[ibody][iv].v.dist(e[triTable[cubeIndex][n+1]]) < 1.0e-3){
		triVertex.y = iv;
		//Add the triangle to this vertex
		bodyVertex[ibody][iv].triangles.push_back(bodyTriangles[ibody].size());
		break;
	      }
	    }

	    if(triVertex.y < 0){
	      triVertex.y = nVertex++;
	      bodyVertex[ibody].push_back(e[triTable[cubeIndex][n+1]]);  //v1
	      //Add the triangle to this vertex
	      bodyVertex[ibody].back().triangles.push_back(bodyTriangles[ibody].size());
	    }
	    

	    //V2
	    for(size_t iv = firstVertexPrev2ZPlane[ibody]; iv < nVertex; ++iv){
	      if(bodyVertex[ibody][iv].v.dist(e[triTable[cubeIndex][n  ]]) < 1.0e-3){
		triVertex.z = iv;
		//Add the triangle to this vertex
		bodyVertex[ibody][iv].triangles.push_back(bodyTriangles[ibody].size());
		break;
	      }
	    }

	    if(triVertex.z < 0){
	      triVertex.z = nVertex++;
	      bodyVertex[ibody].push_back(e[triTable[cubeIndex][n  ]]);  //v2
	      //Add the triangle to this vertex
	      bodyVertex[ibody].back().triangles.push_back(bodyTriangles[ibody].size());
	    }

	    //Save triangle
	    bodyTriangles[ibody].push_back({
		//vertex index
		{static_cast<unsigned long>(triVertex.x),
		static_cast<unsigned long>(triVertex.y),
		static_cast<unsigned long>(triVertex.z)},
		//centroid
		{center},
		//normal
		{0.0,0.0,0.0}
	      });
	  }

	}
      }
    }

    //Store initial number of vertex until two plane before 
    for(unsigned ibody = 0; ibody < nBodies; ++ibody){
      firstVertexPrev2ZPlane[ibody] = firstVertexPrevZPlane[ibody];
    }
    
    //Move planes
    unsigned long* aux = firstPlane;
    firstPlane = lastPlane;
    lastPlane = aux;
    
  }

  printf("\n");
  
  delete [] voxBodies;

  //Post process vertex
  for(unsigned long b = 0; b < nBodies; ++b){

    std::vector<vertex>& bodyV = bodyVertex[b];
    
    //Smooth mesh
    for(unsigned is = 0; is < smoothSteps; ++is){
      for(vertex& v : bodyV){

	//Smooth this vertex
	vector3D<float> midC{0.0f,0.0f,0.0f};
	for(const unsigned it : v.triangles){
	  midC.add(bodyTriangles[b][it].c);
	}
	
	midC = midC/v.triangles.size();
	const vector3D<float> diff = midC-v.v;
	v.v.add(diff*smoothFactor);
      }

      //Recalculate body centers
      for(meshTriangle& t : bodyTriangles[b]){
	t.c = (bodyV[t.vi.x].v +
	       bodyV[t.vi.y].v +
	       bodyV[t.vi.z].v)/3.0f;
      }
    }

    //Compute resulting triangle normals
    for(meshTriangle& t : bodyTriangles[b]){

      t.normal = bodyV[t.vi.y].v - bodyV[t.vi.x].v; //e01
      const vector3D<float> e02 = bodyV[t.vi.z].v   - bodyV[t.vi.x].v;
      t.normal.crossProd(e02);
      t.normal.normalize();
    }    

    //Compute resulting vertex normals
    for(vertex& v : bodyVertex[b]){
      //Calculate the mean normal for this vertex
      vector3D<float> normal = {0.0f,0.0f,0.0f};
      for(const unsigned it : v.triangles){
	normal.add(bodyTriangles[b][it].normal);
      }
      normal.normalize();
      v.normal = normal;
    }
  }
  

  //Print wavefront file
  FILE* fout = nullptr;
  fout = fopen("geo.obj","w");
  if(fout == nullptr){
    printf("Error: Unable to create 'geo.obj' file\n");
    return 1;
  }

  FILE* foutMat = nullptr;
  foutMat = fopen("geoMats.mtl","w");
  if(foutMat == nullptr){
    printf("Error: Unable to create 'geoMats.mtl' file\n");
    return 1;
  }

  //Specify material file
  fprintf(fout,"mtllib geoMats.mtl\n");
  unsigned long nextBodyFirstVertex = 0;
  for(unsigned long b = 0; b < nBodies; ++b){

    unsigned matIndex = b%195;
    
    //Define material for this object
    fprintf(foutMat,"#Material %lu\n",b);
    fprintf(foutMat,"newmtl %lu\n",b);
    fprintf(foutMat,"Ka 1.0 1.0 1.0\n");
    fprintf(foutMat,"Kd %.6f %.6f %.6f\n\n",palette[matIndex*3],palette[matIndex*3+1],palette[matIndex*3+2]);

    fprintf(fout,"usemtl %lu\n",b);
    fprintf(fout,"o Body-%lu\n",b);
    //Vertex coordinates
    for(const vertex& v : bodyVertex[b]){
      fprintf(fout,"v %e %e %e\n", v.v.x, v.v.y, v.v.z);
    }
    //Vertex normals
    for(const vertex& v : bodyVertex[b]){
      fprintf(fout,"vn %e %e %e\n", v.normal.x, v.normal.y, v.normal.z);
    }

    for(const meshTriangle& f : bodyTriangles[b]){
      fprintf(fout,"f %lu//%lu %lu//%lu %lu//%lu\n",
	      f.vi.x+1+nextBodyFirstVertex, f.vi.x+1+nextBodyFirstVertex,
	      f.vi.y+1+nextBodyFirstVertex, f.vi.y+1+nextBodyFirstVertex,
	      f.vi.z+1+nextBodyFirstVertex, f.vi.z+1+nextBodyFirstVertex);
    }
    nextBodyFirstVertex += bodyVertex[b].size();
  }

  fclose(fout);
  fclose(foutMat);
  
  return 0;
}
