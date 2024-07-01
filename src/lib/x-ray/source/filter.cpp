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
//        vicent.gimenez.alventosa@gmail.com
//    
//

#include "filter.hh"

namespace penred{

  namespace xray{

    void createBaseFilter(double dx,
			  double dy,
			  double dz,
			  unsigned nVG,
			  std::ostream& out,
			  const unsigned matIndex,
			  const std::string& filterName,
			  const std::string& parentName,
			  const bool numObjects,
			  const vector3D<double> center){
      
      //Creates a deformable mesh filter with dimensions "dx" x "dy" and
      //"nVG" vertex groups uniformly distributed on the x axis.
      //
      //
      //  Vertex are constructed as x slices from -x to +x with the following order: 
      //
      // 3*---------------------*2  n Z
      //  |                     |   |
      //  |                     |   |
      // 0*---------------------*1  |------> Y
      //
      // i.e. counterclock in the y-z plane. For example, for 8 vertex groups on
      // the "x" axis viewed from above in the Z direction, the vertex numering are
      // as follows:
      //                             
      // 2  6 10 14 18 22 26 30   n Y
      // *--*--*--*--*--*--*--*   | 
      // |                    |   |
      // *--*--*--*--*--*--*--*   |-----> X
      // 3  7 11 15 19 23 27 31
      //
      //
      // with a total of nVG * 4 = 32 vertex

      //Check dimensions
      if(dx <= 0.0)
	dx = 1.0;
      if(dy <= 0.0)
	dy = 1.0;
      if(dz <= 0.0)
	dz = 1.0;
      
      //Ensure to have, at least, two vertex groups in the x axis
      nVG = std::max(nVG,2u);
      
      //Calculate vertex distance in X axis
      const double dvg = dx/static_cast<double>(nVG-1);

      //Calculate number of vertex
      const unsigned nVertex = nVG*4;

      //Calculate number of triangles
      const unsigned nFaces = 8*(nVG-1) + 4; //+4 for lateral faces on "x" axis

      //Create the resulting string
      std::string result;
      
      //Check if the number of objects must be printed
      if(numObjects){
	out << "# Number of objects:\n 1" << std::endl;
      }

      // ** Print object name and header
      out << "# Object: " << filterName << std::endl;
      out << "#MAT      #NFACES     #NVERTEX     #NAME"
	"        #PARENT NAME    #N VERTEX GROUPS"
	  << std::endl;;
      
      //Material index
      out << " " << std::to_string(matIndex);

      //Number of faces
      out << "   " << std::to_string(nFaces);
      
      //Number of vertex
      out << "   " << std::to_string(nVertex);

      //Filter name
      out << "   " << filterName;

      //Parent name
      out << "   " << parentName;

      //Vertex groups. Add 2 for "y" axis limits and 2 for "z" limits
      // Also, "x" groups are split in top and bottom vertex
      out << "   " <<  2*nVG + 4 << std::endl;

      // ** Print vertex groups
      out << "# VERTEX GROUPS" << std::endl;

      //+y group
      out << "#NAME  #NVERTEX"<< std::endl;
      out << " +y " << 2*nVG << std::endl;
      for(size_t ivg = 0; ivg < nVG; ++ivg){
	unsigned iStart = ivg*4; // 4 vertex per group
	out << iStart + 1 << std::endl;
	out << iStart + 2 << std::endl;
      }

      //-y group
      out << "#NAME  #NVERTEX"<< std::endl;
      out << " -y " << 2*nVG << std::endl;
      for(size_t ivg = 0; ivg < nVG; ++ivg){
	unsigned iStart = ivg*4; // 4 vertex per group
	out << iStart     << std::endl;
	out << iStart + 3 << std::endl;
      }

      //+z group
      out << "#NAME  #NVERTEX"<< std::endl;
      out << " +z " << 2*nVG << std::endl;
      for(size_t ivg = 0; ivg < nVG; ++ivg){
	unsigned iStart = ivg*4; // 4 vertex per group
	out << iStart + 2 << std::endl;
	out << iStart + 3 << std::endl;
      }

      //-z group
      out << "#NAME  #NVERTEX"<< std::endl;
      out << " -z " << 2*nVG << std::endl;
      for(size_t ivg = 0; ivg < nVG; ++ivg){
	unsigned iStart = ivg*4; // 4 vertex per group
	out << iStart     << std::endl;
	out << iStart + 1 << std::endl;
      }      

      //Create X groups with positive Z
      for(size_t ivg = 0; ivg < nVG; ++ivg){
	out << "#NAME  #NVERTEX"<< std::endl;
	out << " +x-" << ivg << " 2" << std::endl;
	unsigned iStart = ivg*4; // 4 vertex per group
	out << iStart + 2 << std::endl;
	out << iStart + 3 << std::endl;
      }

      //Create X groups with negative Z
      for(size_t ivg = 0; ivg < nVG; ++ivg){
	out << "#NAME  #NVERTEX"<< std::endl;
	out << " -x-" << ivg << " 2" << std::endl;
	unsigned iStart = ivg*4; // 4 vertex per group
	out << iStart     << std::endl;
	out << iStart + 1 << std::endl;
      }      

      // ** Print vertex positions
      out << "# VERTEX LIST" << std::endl;
      out << "# Index (X Y Z)" << std::endl;

      const double originX = center.x-dx/2.0;
      const double originY = center.y-dy/2.0;
      const double originZ = center.z-dz/2.0;
      const double limitY = center.y+dy/2.0;
      const double limitZ = center.z+dz/2.0;
      for(size_t ivg = 0; ivg < nVG; ++ivg){
	double startIndex = ivg*4;
	double xShift = static_cast<double>(ivg)*dvg;
	double xPos = originX + xShift;

	//Print the four vertex in this group

	// 0
	out << startIndex
	    << " " << xPos << " " << originY << " " << originZ << std::endl;

	// 1
	out << startIndex + 1
	    << " " << xPos << " " << limitY  << " " << originZ << std::endl;

	// 2
	out << startIndex + 2
	    << " " << xPos << " " << limitY  << " " << limitZ  << std::endl;

	// 3
	out << startIndex + 3
	    << " " << xPos << " " << originY << " " << limitZ  << std::endl;	
      }

      // ** Construct triangles
      // take into account the triangle orientation to avoid erroneous face culling
      out << "# FACES(triangles)" << std::endl;
      for(size_t ivg = 0; ivg < nVG-1; ++ivg){
	unsigned long startIndex = ivg*4;

	//Z bot triangles
	//
	//  4                     5   n X
	//  *---------------------*   |
	//  |                     |   |
	//  |                     |   |
	//  *---------------------*   |--------> Y
	//  0                     1

	// 0-1-5
	out << startIndex     << " "
	    << startIndex + 1 << " "
	    << startIndex + 5 << std::endl;
	// 0-5-4
	out << startIndex     << " "
	    << startIndex + 5 << " "
	    << startIndex + 4 << std::endl;


	//Z top triangles
	//
	//  7                     6   n X
	//  *---------------------*   |
	//  |                     |   |
	//  |                     |   |
	//  *---------------------*   |--------> Y
	//  3                     2

	// 3-6-2
	out << startIndex + 3 << " "
	    << startIndex + 6 << " "
	    << startIndex + 2 << std::endl;
	// 3-7-6
	out << startIndex + 3 << " "
	    << startIndex + 7 << " "
	    << startIndex + 6 << std::endl;

	//Y bot triangles
	//
	//  3      7   n Z
	//  *------*   |
	//  |      |   |
	//  |      |   |
	//  *------*   |--------> X
	//  0      4

	// 0-4-7
	out << startIndex     << " "
	    << startIndex + 4 << " "
	    << startIndex + 7 << std::endl;
	// 0-7-3
	out << startIndex     << " "
	    << startIndex + 7 << " "
	    << startIndex + 3 << std::endl;

	
	//Y top triangles
	//
	//  2      6   n Z
	//  *------*   |
	//  |      |   |
	//  |      |   |
	//  *------*   |--------> X
	//  1      5

	// 1-6-5
	out << startIndex + 1 << " "
	    << startIndex + 6 << " "
	    << startIndex + 5 << std::endl;
	// 1-2-6
	out << startIndex + 1 << " "
	    << startIndex + 2 << " "
	    << startIndex + 6 << std::endl;
	
      }

      // ** Create faces at X limits

      //
      // 3*---------------------*2  n Z
      //  |                     |   |
      //  |                     |   |
      // 0*---------------------*1  |------> Y

      //X bot triangles
      // 0-2-1
      out << " 0 2 1" << std::endl;
      // 0-3-2
      out << " 0 3 2" << std::endl;

      //X top triangles
      // 0-1-2
      unsigned long startIndex = (nVG-1)*4;
      out << startIndex     << " "
	  << startIndex + 1 << " "
	  << startIndex + 2 << std::endl;
      // 0-2-3
      out << startIndex     << " "
	  << startIndex + 2 << " "
	  << startIndex + 3 << std::endl;
      
    }

    void createTopFaceIrregularFilter(double dx,
				      double dy,
				      std::vector<double>& dz,
				      std::ostream& out,
				      const unsigned matIndex,
				      const std::string& filterName,
				      const std::string& parentName,
				      const bool numObjects,
				      const vector3D<double> center){
      
      //Creates a deformable mesh filter with dimensions "dx" x "dy" and
      //width and number of vertex groups depending on the "dz" variable.
      //Vertex groups are uniformly distributed on the x axis.
      //
      //
      //  Vertex are constructed as x slices from -x to +x with the following order: 
      //
      // 3*---------------------*2  n Z
      //  |                     |   |
      //  |                     |   |
      // 0*---------------------*1  |------> Y
      //
      // i.e. counterclock in the y-z plane. For example, for 8 vertex groups on
      // the "x" axis viewed from above in the Z direction, the vertex numering are
      // as follows:
      //                             
      // 2  6 10 14 18 22 26 30   n Y
      // *--*--*--*--*--*--*--*   | 
      // |                    |   |
      // *--*--*--*--*--*--*--*   |-----> X
      // 3  7 11 15 19 23 27 31
      //
      //
      // with a total of nVG * 4 = 32 vertex

      //Check dimensions
      if(dx <= 0.0)
	dx = 1.0;
      if(dy <= 0.0)
	dy = 1.0;

      double maxDz = -1.0;
      for(double& z : dz){
	if(z <= 0.0)
	  z = 0.01;
	if(maxDz < z){
	  maxDz = z;
	}
      }

      if(dz.size() < 2){
	createBaseFilter(dx, dy, dz[0], 0, out,
			 matIndex, filterName, parentName,
			 numObjects, center);
	return;
      }
      
      //Ensure to have, at least, two vertex groups in the x axis
      size_t nVG = dz.size();
      
      //Calculate vertex distance in X axis
      const double dvg = dx/static_cast<double>(nVG-1);

      //Calculate number of vertex
      const unsigned nVertex = nVG*4;

      //Calculate number of triangles
      const unsigned nFaces = 8*(nVG-1) + 4; //+4 for lateral faces on "x" axis

      //Create the resulting string
      std::string result;
      
      //Check if the number of objects must be printed
      if(numObjects){
	out << "# Number of objects:\n 1" << std::endl;
      }

      // ** Print object name and header
      out << "# Object: " << filterName << std::endl;
      out << "#MAT      #NFACES     #NVERTEX     #NAME"
	"        #PARENT NAME    #N VERTEX GROUPS"
	  << std::endl;;
      
      //Material index
      out << " " << std::to_string(matIndex);

      //Number of faces
      out << "   " << std::to_string(nFaces);
      
      //Number of vertex
      out << "   " << std::to_string(nVertex);

      //Filter name
      out << "   " << filterName;

      //Parent name
      out << "   " << parentName;

      //Vertex groups. Add 2 for "y" axis limits and 2 for "z" limits
      // Also, "x" groups are split in top and bottom vertex
      out << "   " <<  2*nVG + 4 << std::endl;

      // ** Print vertex groups
      out << "# VERTEX GROUPS" << std::endl;

      //+y group
      out << "#NAME  #NVERTEX"<< std::endl;
      out << " +y " << 2*nVG << std::endl;
      for(size_t ivg = 0; ivg < nVG; ++ivg){
	unsigned iStart = ivg*4; // 4 vertex per group
	out << iStart + 1 << std::endl;
	out << iStart + 2 << std::endl;
      }

      //-y group
      out << "#NAME  #NVERTEX"<< std::endl;
      out << " -y " << 2*nVG << std::endl;
      for(size_t ivg = 0; ivg < nVG; ++ivg){
	unsigned iStart = ivg*4; // 4 vertex per group
	out << iStart     << std::endl;
	out << iStart + 3 << std::endl;
      }

      //+z group
      out << "#NAME  #NVERTEX"<< std::endl;
      out << " +z " << 2*nVG << std::endl;
      for(size_t ivg = 0; ivg < nVG; ++ivg){
	unsigned iStart = ivg*4; // 4 vertex per group
	out << iStart + 2 << std::endl;
	out << iStart + 3 << std::endl;
      }

      //-z group
      out << "#NAME  #NVERTEX"<< std::endl;
      out << " -z " << 2*nVG << std::endl;
      for(size_t ivg = 0; ivg < nVG; ++ivg){
	unsigned iStart = ivg*4; // 4 vertex per group
	out << iStart     << std::endl;
	out << iStart + 1 << std::endl;
      }      

      //Create X groups with positive Z
      for(size_t ivg = 0; ivg < nVG; ++ivg){
	out << "#NAME  #NVERTEX"<< std::endl;
	out << " +x-" << ivg << " 2" << std::endl;
	unsigned iStart = ivg*4; // 4 vertex per group
	out << iStart + 2 << std::endl;
	out << iStart + 3 << std::endl;
      }

      //Create X groups with negative Z
      for(size_t ivg = 0; ivg < nVG; ++ivg){
	out << "#NAME  #NVERTEX"<< std::endl;
	out << " -x-" << ivg << " 2" << std::endl;
	unsigned iStart = ivg*4; // 4 vertex per group
	out << iStart     << std::endl;
	out << iStart + 1 << std::endl;
      }      

      // ** Print vertex positions
      out << "# VERTEX LIST" << std::endl;
      out << "# Index (X Y Z)" << std::endl;

      const double originX = center.x-dx/2.0;
      const double originY = center.y-dy/2.0;
      const double originZ = center.z-maxDz/2.0;
      const double limitY = center.y+dy/2.0;
      for(size_t ivg = 0; ivg < nVG; ++ivg){
	double startIndex = ivg*4;
	double xShift = static_cast<double>(ivg)*dvg;
	double xPos = originX + xShift;
	const double finalZ = originZ + dz[ivg];

	//Print the four vertex in this group

	// 0
	out << startIndex
	    << " " << xPos << " " << originY << " " << originZ << std::endl;

	// 1
	out << startIndex + 1
	    << " " << xPos << " " << limitY  << " " << originZ << std::endl;

	// 2
	out << startIndex + 2
	    << " " << xPos << " " << limitY  << " " << finalZ  << std::endl;

	// 3
	out << startIndex + 3
	    << " " << xPos << " " << originY << " " << finalZ  << std::endl;	
      }

      // ** Construct triangles
      // take into account the triangle orientation to avoid erroneous face culling
      out << "# FACES(triangles)" << std::endl;
      for(size_t ivg = 0; ivg < nVG-1; ++ivg){
	unsigned long startIndex = ivg*4;

	//Z bot triangles
	//
	//  4                     5   n X
	//  *---------------------*   |
	//  |                     |   |
	//  |                     |   |
	//  *---------------------*   |--------> Y
	//  0                     1

	// 0-1-5
	out << startIndex     << " "
	    << startIndex + 1 << " "
	    << startIndex + 5 << std::endl;
	// 0-5-4
	out << startIndex     << " "
	    << startIndex + 5 << " "
	    << startIndex + 4 << std::endl;


	//Z top triangles
	//
	//  7                     6   n X
	//  *---------------------*   |
	//  |                     |   |
	//  |                     |   |
	//  *---------------------*   |--------> Y
	//  3                     2

	// 3-6-2
	out << startIndex + 3 << " "
	    << startIndex + 6 << " "
	    << startIndex + 2 << std::endl;
	// 3-7-6
	out << startIndex + 3 << " "
	    << startIndex + 7 << " "
	    << startIndex + 6 << std::endl;

	//Y bot triangles
	//
	//  3      7   n Z
	//  *------*   |
	//  |      |   |
	//  |      |   |
	//  *------*   |--------> X
	//  0      4

	// 0-4-7
	out << startIndex     << " "
	    << startIndex + 4 << " "
	    << startIndex + 7 << std::endl;
	// 0-7-3
	out << startIndex     << " "
	    << startIndex + 7 << " "
	    << startIndex + 3 << std::endl;

	
	//Y top triangles
	//
	//  2      6   n Z
	//  *------*   |
	//  |      |   |
	//  |      |   |
	//  *------*   |--------> X
	//  1      5

	// 1-6-5
	out << startIndex + 1 << " "
	    << startIndex + 6 << " "
	    << startIndex + 5 << std::endl;
	// 1-2-6
	out << startIndex + 1 << " "
	    << startIndex + 2 << " "
	    << startIndex + 6 << std::endl;
	
      }

      // ** Create faces at X limits

      //
      // 3*---------------------*2  n Z
      //  |                     |   |
      //  |                     |   |
      // 0*---------------------*1  |------> Y

      //X bot triangles
      // 0-2-1
      out << " 0 2 1" << std::endl;
      // 0-3-2
      out << " 0 3 2" << std::endl;

      //X top triangles
      // 0-1-2
      unsigned long startIndex = (nVG-1)*4;
      out << startIndex     << " "
	  << startIndex + 1 << " "
	  << startIndex + 2 << std::endl;
      // 0-2-3
      out << startIndex     << " "
	  << startIndex + 2 << " "
	  << startIndex + 3 << std::endl;
      
    }

  } // namespace xray
} // namespace penred
