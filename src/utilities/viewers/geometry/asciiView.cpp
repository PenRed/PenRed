//
//
//    Copyright (C) 2021 Universitat de València - UV
//    Copyright (C) 2021 Universitat Politècnica de València - UPV
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
//        vicent.gimenez.alventosa@gmail.com  (Vicent Giménez Alventosa)
//        vicente.gimenez@uv.es (Vicent Giménez Gómez)
//        sanolgi@upvnet.upv.es (Sandra Olver Gil)
//

#include <array>
#include "pen_geoView.hh"

void printRender(std::string filename, unsigned char* render, unsigned nx, unsigned ny);

int main(int argc, char** argv){

  if(argc < 2){
    printf("usage: %s configure-file\n",argv[0]);
    return 1;
  }

  //Parse configuration file
  pen_parserSection config;
  std::string errorLine;
  unsigned long errorLineNum;
  int err = parseFile(argv[1],config,errorLine,errorLineNum);
  
  if(err != INTDATA_SUCCESS){
    printf("Error parsing configuration.\n");
    printf("Error code: %d\n",err);
    printf("Error message: %s\n",pen_parserError(err));
    printf("Error located at line %lu, at text: %s\n",
	   errorLineNum,errorLine.c_str());
    return -1;
  }
  
  //Get geometry section
  pen_parserSection geometrySection;
  err = config.readSubsection("geometry",geometrySection);
  if(err != INTDATA_SUCCESS){
    printf("Unable to read 'geometry' section.\n");
    printf("Error code: %d\n",err);
    printf("Error message: %s\n",pen_parserError(err));
    return -2;
  }

  pen_geoView viewer;
  err = viewer.init(geometrySection,3);
  if(err != 0){
    printf("Error: Unable to initialize viewer");
    return -3;
  }

  const unsigned nxmax = 512;
  const unsigned nymax = 512;
  const unsigned nxymax = nxmax*nymax;
  std::array<unsigned char, nxymax> renderMat;
  std::array<unsigned int, nxymax> renderBody;
  std::array<float, nxymax> distances;

  //2D renders
  //Extract 2D section
  pen_parserSection section2D;
  if(config.readSubsection("2D",section2D) != INTDATA_SUCCESS){
    printf("No 2D planes specified.\n");
  }
  else{
    printf("Printing specified 2D planes.\n");

    // Z planes
    //************
    pen_parserSection section2Dz;
    if(section2D.readSubsection("z",section2Dz) != INTDATA_SUCCESS){
      printf("No 2D z planes specified.\n");
    }
    else{
      std::vector<std::string> names;
      section2Dz.ls(names);

      printf("Number of specified z planes: %u\n",static_cast<unsigned>(names.size()));
      
      for(size_t i = 0; i < names.size(); ++i){

	//Read plane section
	pen_parserSection planeSec;
	if(section2Dz.read(names[i],planeSec) != INTDATA_SUCCESS)
	  continue;
	
	pen_parserArray pos;
	if(planeSec.read("pos",pos) != INTDATA_SUCCESS){
	  printf("Unable to read position ('pos') of z plane '%s'\n",names[i].c_str());
	}
	else{
	  double x,y,z;
	  x = y = z = 0.0;
	  pos[0].read(x);
	  pos[1].read(y);
	  pos[2].read(z);
	  
	  double dsPixel;
	  if(planeSec.read("pixel-size",dsPixel) != INTDATA_SUCCESS){
	    printf("Unable to read pixel size ('pixel-size') of z plane '%s'\n",names[i].c_str());
	  }
	  else{

	    //Read optional values nx and ny
	    int nx = nxmax;
	    int ny = nymax;
	    planeSec.read("nx",nx);
	    planeSec.read("ny",ny);

	    if(nx <= 0 || ny <= 0){
	      printf("Invalid number of pixels.\n"
		     "    Pixels in X axis: %d\n"
		     "    Pixels in Y axis: %d\n",nx,ny);
	    }
	    else{
	      viewer.renderZ(renderMat.data(),renderBody.data(),
			     x,y,z,
			     dsPixel,dsPixel,
			     static_cast<unsigned>(nx),static_cast<unsigned>(ny));


	      std::string filename = names[i] + std::string("_z_") + std::to_string(dsPixel) +
		std::string("_") + std::to_string(nx) + std::string("x") + std::to_string(ny)
		+ std::string(".dat");
	      printRender(filename, renderMat.data(), static_cast<unsigned>(nx),static_cast<unsigned>(ny));
	    }
	  }
	}
      }
      
    }

    // Y planes
    //************
    pen_parserSection section2Dy;
    if(section2D.readSubsection("y",section2Dy) != INTDATA_SUCCESS){
      printf("No 2D y planes specified.\n");
    }
    else{
      std::vector<std::string> names;
      section2Dy.ls(names);

      printf("Number of specified y planes: %u\n",static_cast<unsigned>(names.size()));
      
      for(size_t i = 0; i < names.size(); ++i){

	//Read plane section
	pen_parserSection planeSec;
	if(section2Dy.read(names[i],planeSec) != INTDATA_SUCCESS)
	  continue;
	
	pen_parserArray pos;
	if(planeSec.read("pos",pos) != INTDATA_SUCCESS){
	  printf("Unable to read position ('pos') of y plane '%s'\n",names[i].c_str());
	}
	else{
	  double x,y,z;
	  x = y = z = 0.0;
	  pos[0].read(x);
	  pos[1].read(y);
	  pos[2].read(z);
	  
	  double dsPixel;
	  if(planeSec.read("pixel-size",dsPixel) != INTDATA_SUCCESS){
	    printf("Unable to read pixel size ('pixel-size') of y plane '%s'\n",names[i].c_str());
	  }
	  else{

	    //Read optional values nx and ny
	    int nx = nxmax;
	    int nz = nymax;
	    planeSec.read("nx",nx);
	    planeSec.read("nz",nz);

	    if(nx <= 0 || nz <= 0){
	      printf("Invalid number of pixels.\n"
		     "    Pixels in X axis: %d\n"
		     "    Pixels in Z axis: %d\n",nx,nz);
	    }
	    else{
	      viewer.renderY(renderMat.data(),renderBody.data(),
			     x,y,z,
			     dsPixel,dsPixel,
			     static_cast<unsigned>(nx),static_cast<unsigned>(nz));


	      std::string filename = names[i] + std::string("_y_") + std::to_string(dsPixel) +
		std::string("_") + std::to_string(nx) + std::string("x") + std::to_string(nz)
		+ std::string(".dat");
	      printRender(filename, renderMat.data(), static_cast<unsigned>(nx),static_cast<unsigned>(nz));
	    }
	  }
	}
      }
      
    }

    // X planes
    //************
    pen_parserSection section2Dx;
    if(section2D.readSubsection("x",section2Dx) != INTDATA_SUCCESS){
      printf("No 2D x planes specified.\n");
    }
    else{
      std::vector<std::string> names;
      section2Dx.ls(names);

      printf("Number of specified x planes: %u\n",static_cast<unsigned>(names.size()));
      
      for(size_t i = 0; i < names.size(); ++i){

	//Read plane section
	pen_parserSection planeSec;
	if(section2Dx.read(names[i],planeSec) != INTDATA_SUCCESS)
	  continue;
	
	pen_parserArray pos;
	if(planeSec.read("pos",pos) != INTDATA_SUCCESS){
	  printf("Unable to read position ('pos') of x plane '%s'\n",names[i].c_str());
	}
	else{
	  double x,y,z;
	  x = y = z = 0.0;
	  pos[0].read(x);
	  pos[1].read(y);
	  pos[2].read(z);
	  
	  double dsPixel;
	  if(planeSec.read("pixel-size",dsPixel) != INTDATA_SUCCESS){
	    printf("Unable to read pixel size ('pixel-size') of x plane '%s'\n",names[i].c_str());
	  }
	  else{

	    //Read optional values nx and ny
	    int ny = nxmax;
	    int nz = nymax;
	    planeSec.read("ny",ny);
	    planeSec.read("nz",nz);

	    if(ny <= 0 || nz <= 0){
	      printf("Invalid number of pixels.\n"
		     "    Pixels in Y axis: %d\n"
		     "    Pixels in Z axis: %d\n",ny,nz);
	    }
	    else{
	      viewer.renderX(renderMat.data(),renderBody.data(),
			     x,y,z,
			     dsPixel,dsPixel,
			     static_cast<unsigned>(ny),static_cast<unsigned>(nz));


	      std::string filename = names[i] + std::string("_x_") + std::to_string(dsPixel) +
		std::string("_") + std::to_string(ny) + std::string("x") + std::to_string(nz)
		+ std::string(".dat");
	      printRender(filename, renderMat.data(), static_cast<unsigned>(ny),static_cast<unsigned>(nz));
	    }
	  }
	}
      }
      
    }
    
  }

  //3D renders
  //Extract 3D section
  pen_parserSection section3D;
  if(config.readSubsection("3D",section3D) != INTDATA_SUCCESS){
    printf("No 3D positions specified.\n");
  }
  else{
    printf("Printing specified 3D positions.\n");

    std::vector<std::string> names;
    section3D.ls(names);

    for(size_t i = 0; i < names.size(); ++i){

	//Read plane section
	pen_parserSection planeSec;
	if(section3D.read(names[i],planeSec) != INTDATA_SUCCESS)
	  continue;
	
	pen_parserArray pos;
	if(planeSec.read("pos",pos) != INTDATA_SUCCESS){
	  printf("Unable to read position ('pos') of 3D render '%s'\n",names[i].c_str());
	}
	else{
	  double x,y,z;
	  x = y = z = 0.0;
	  pos[0].read(x);
	  pos[1].read(y);
	  pos[2].read(z);

	  pen_parserArray dir;
	  if(planeSec.read("dir",dir) != INTDATA_SUCCESS){
	    printf("Unable to read direction ('dir') of 3D render '%s'\n",names[i].c_str());
	  }
	  else{

	    double u,v,w;
	    u = v = 0.0;
	    w = 1.0;
	    dir[0].read(u);
	    dir[1].read(v);
	    dir[2].read(w);
	    
	    double dsPixel;
	    if(planeSec.read("pixel-size",dsPixel) != INTDATA_SUCCESS){
	      printf("Unable to read pixel size ('pixel-size') of x plane '%s'\n",names[i].c_str());
	    }
	    else{

	      //Read optional values nx and ny
	      int nx = nxmax;
	      int ny = nymax;
	      planeSec.read("nx",nx);
	      planeSec.read("ny",ny);

	      if(nx <= 0 || ny <= 0){
		printf("Invalid number of pixels.\n"
		       "    Pixels in X axis: %d\n"
		       "    Pixels in Y axis: %d\n",nx,ny);
	      }
	      else{

		double roll;
		planeSec.read("roll",roll);

		float phi = 0.0;
		float minD,maxD;
		viewer.render3D(renderMat.data(),renderBody.data(),
				x,y,z,
				u,v,w,
				roll, phi, // auxiliar phi
				distances.data(),minD,maxD);
		//dsPixel,dsPixel,
		//static_cast<unsigned>(nx),static_cast<unsigned>(ny));


		std::string filename = names[i] + std::string("_3D_") + std::to_string(dsPixel) +
		  std::string("_") + std::to_string(nx) + std::string("x") + std::to_string(ny)
		  + std::string(".dat");
		printRender(filename, renderMat.data(), static_cast<unsigned>(nx),static_cast<unsigned>(ny));
	      }
	    }
	  }
	}
      
    }
  }
  
  return 0;
  
}

void printRender(std::string filename, unsigned char* render, unsigned nx, unsigned ny){

  FILE* fout = nullptr;
  fout = fopen(filename.c_str(),"w");
  if(fout == nullptr)
    return;
  
  for(unsigned j = 0; j < ny; ++j){
    for(unsigned i = 0; i < nx; ++i){
      fprintf(fout,"%u ", static_cast<unsigned>(render[j*nx+i]));
    }
    fprintf(fout,"\n");
  }
  
}

