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

#include "collimator.hh"

namespace penred{

  namespace xray{

    bool collimateInVoid(pen_particleState& state,
			 const double zUp,
			 const double zDown,
			 const double dxUp, const double dyUp,
			 const double dxDown, const double dyDown){

      //Reject particles going up
      if(state.W >= 0.0){
	return false;
      }

      //Check if the particle is under the collimator
      if(state.Z <= zDown){
	//Do not change this particle state
	return true;
      }
	
      //Check if the particle is
      if(state.Z > zUp){
	//Particle is above the collimator
	  
	//Get the distance until the collimato
	double ds2Col = (state.Z-zUp)/state.W;

	//Check if the particle reaches the enter limits
	double xIn = state.X + state.U*ds2Col;
	double yIn = state.Y + state.V*ds2Col;

	if(std::fabs(xIn) >= dxUp/2.0 ||
	   std::fabs(yIn) >= dyUp/2.0){
	  return false;
	}
      }
	
      //At this point, the particle is inside, or can enter, the collimator region
      //Check if the limits are reached at the collimator end
      double ds2Out = (state.Z-zDown)/state.W;
      double xOut = state.X + state.U*ds2Out;
      double yOut = state.Y + state.V*ds2Out;

      if(std::fabs(xOut) >= dxDown/2.0 ||
	 std::fabs(yOut) >= dyDown/2.0){
	return false;
      }
	
      //Particle cross the collimator.
      //Move it until the collimator end
      state.X = xOut;
      state.Y = yOut;
      state.Z = zDown;

      return true;
    }

    void collimateInVoid(std::vector<detectedPart>& beam,
			 const double zUp,
			 const double zDown,
			 const double dxUp, const double dyUp,
			 const double dxDown, const double dyDown){

      //Iterate over particles
      size_t i = 0;
      size_t nPart = beam.size();
      //Accumulate history increment of skiped particles
      unsigned long long dh = 0; 
      while(i < nPart){
	
	//Get particle state
	pen_particleState& state = beam[i].state;

	//Reject particles going up
	if(state.W >= 0.0){
	  //Accumulate history increment
	  dh += beam[i].dh;
	  beam[i] = beam[--nPart];
	  continue;
	}

	//Check if the particle is under the collimator
	if(state.Z <= zDown){
	  //This particle pass the collimator.
	  //Add history increment and restart it
	  beam[i].dh += dh;
	  dh = 0;
	  ++i;
	  continue;
	}
	
	//Check if the particle is
	if(state.Z > zUp){
	  //Particle is above the collimator
	  
	  //Get the distance until the collimato
	  double ds2Col = (state.Z-zUp)/state.W;

	  //Check if the particle reaches the enter limits
	  double xIn = state.X + state.U*ds2Col;
	  double yIn = state.Y + state.V*ds2Col;

	  if(std::fabs(xIn) >= dxUp/2.0 ||
	     std::fabs(yIn) >= dyUp/2.0){
	    //Accumulate history increment
	    dh += beam[i].dh;    
	    //Reject the particle	    
	    beam[i] = beam[--nPart];
	    continue;
	  }
	}
	
	//At this point, the particle is inside, or can enter, the collimator region
	//Check if the limits are reached at the collimator end
	double ds2Out = (state.Z-zDown)/state.W;
	double xOut = state.X + state.U*ds2Out;
	double yOut = state.Y + state.V*ds2Out;

	if(std::fabs(xOut) >= dxDown/2.0 ||
	   std::fabs(yOut) >= dyDown/2.0){
	  //Accumulate history increment
	  dh += beam[i].dh;	  
	  //Reject the particle
	  beam[i] = beam[--nPart];
	  continue;
	}
	
	//Particle cross the collimator.
	//Move it until the collimator end
	state.X = xOut;
	state.Y = yOut;
	state.Z = zDown;

	//Add history increment and restart it
	beam[i].dh += dh;
	dh = 0;
	
	++i;
      }

      //Resize final beam vector
      beam.resize(nPart);
    }

    
    void createBaseCollimator(double dx, double dxInTop, double dxInBot,
			      double dy, double dyInTop, double dyInBot,
			      double dz,
			      std::ofstream& out,
			      const unsigned matIndex,
			      const std::string& collimatorName,
			      const std::string& parentName,
			      const bool numObjects){
      
      //Create a rectangular shaped collimator as a deformable mesh with initial
      //dimensions "dx" x "dy" x "dz". The collimator is constructed using 4 external
      //vertex for the external bottom part (Z min), 4 for the internal, and the same
      //for the top part (Z top). The vertex indexes are assigned in the range [0,16)
      //as is shown in the figures below, and are symetrics for the bot and top parts.
      //
      //       Bottom view
      //
      //  3 *----------------* 2
      //    |\              /|                            -y view
      //    | \            / |
      //    | 7*----------*6 |
      //    |  |          |  |      n Y            8                      9      n Z
      //    |  |   Empty  |  |  dy  |              *----------------------*      |
      //    |  |          |  |      |              |                      | dz   |
      //    | 4*----------*5 |      |              *----------------------*      |
      //    | /            \ |      |--------> X   0           dx         1      |--------> X
      //    |/              \|
      //  0 *----------------* 1
      //            dx
      //
      //
      //
      //        Top view
      //
      // 11 *----------------* 10
      //    |\              /|
      //    | \            / |
      //    |15*----------*14|
      //    |  |          |  |      n Y            
      //    |  |   Empty  |  |  dy  |              
      //    |  |          |  |      |              
      //    |12*----------*13|      |              
      //    | /            \ |      |--------> X  
      //    |/              \|
      //  8 *----------------* 9
      //            dx

      
      //Check dimensions
      if(dx <= 0.0)
	dx = 1.0;
      if(dy <= 0.0)
	dy = 1.0;
      if(dz <= 0.0)
	dz = 1.0;

      if(dxInTop >= dx)
	dxInTop = dx/2.0;
      if(dyInTop >= dy)
	dyInTop = dy/2.0;

      if(dxInBot >= dx)
	dxInBot = dx/2.0;
      if(dyInBot >= dy)
	dyInBot = dy/2.0;
      
      //Set number of vertex
      const unsigned nVertex = 16;

      //Set number of triangles: 8*2 for bot and top limits,
      //4*2 for outer lateral faces and 4*2 for inner lateral faces
      constexpr unsigned nFaces = 8*2 + 4*2 + 2*4;

      //Check if the number of objects must be printed
      if(numObjects){
	out << "# Number of objects:\n 1" << std::endl;
      }

      // ** Print object name and header
      out << "# Object: " << collimatorName << std::endl;
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
      out << "   " << collimatorName;

      //Parent name
      out << "   " << parentName;

      //Number of vertex groups: 4 for lateral vertex, for inner and outer vertex,
      //and 2 for top and for bottom vertex
      out << " " << 4 * 2 + 2 << std::endl;

      // ** Print vertex groups
      out << "# VERTEX GROUPS" << std::endl;

      // External vertex groups
      
      //-y external group
      out << "#NAME  #NVERTEX"<< std::endl;
      out << " -y-out 4" << std::endl;
      out << " 0\n 1\n 8\n 9" << std::endl;

      //+x external group
      out << "#NAME  #NVERTEX"<< std::endl;
      out << " +x-out 4" << std::endl;
      out << " 1\n 2\n 9\n 10" << std::endl;

      //+y external group
      out << "#NAME  #NVERTEX"<< std::endl;
      out << " +y-out 4" << std::endl;
      out << " 2\n 3\n 10\n 11" << std::endl;
      
      //-x external group
      out << "#NAME  #NVERTEX"<< std::endl;
      out << " -x-out 4" << std::endl;
      out << " 3\n 0\n 11\n 8" << std::endl;

      // Internal vertex groups
      
      //-y external group
      out << "#NAME  #NVERTEX"<< std::endl;
      out << " -y-in 4" << std::endl;
      out << " 4\n 5\n 12\n 13" << std::endl;

      //+x external group
      out << "#NAME  #NVERTEX"<< std::endl;
      out << " +x-in 4" << std::endl;
      out << " 5\n 6\n 13\n 14" << std::endl;

      //+y external group
      out << "#NAME  #NVERTEX"<< std::endl;
      out << " +y-in 4" << std::endl;
      out << " 6\n 7\n 14\n 15" << std::endl;
      
      //-x external group
      out << "#NAME  #NVERTEX"<< std::endl;
      out << " -x-in 4" <<  std::endl;
      out << " 7\n 4\n 15\n 12" << std::endl;

      // Bot vertex group
      out << "#NAME  #NVERTEX"<< std::endl;
      out << " bot 8" << std::endl;
      out << " 0\n 1\n 2\n 3\n 4\n 5\n 6\n 7" << std::endl;

      // Top vertex group
      out << "#NAME  #NVERTEX"<< std::endl;
      out << " top 8" << std::endl;
      out << " 8\n 9\n 10\n 11\n 12\n 13\n 14\n 15" << std::endl;


      // ** Print vertex positions
      out << "# VERTEX LIST" << std::endl;
      out << "# Index (X Y Z)" << std::endl;

      const double originX = -dx/2.0;
      const double originY = -dy/2.0;

      const double originInXTop = -dxInTop/2.0;
      const double originInYTop = -dyInTop/2.0;

      const double originInXBot = -dxInBot/2.0;
      const double originInYBot = -dyInBot/2.0;
      
      const double limitX = dx/2.0;
      const double limitY = dy/2.0;

      const double limitInXTop = dxInTop/2.0;
      const double limitInYTop = dyInTop/2.0;

      const double limitInXBot = dxInBot/2.0;
      const double limitInYBot = dyInBot/2.0;
      
      const double originZ = -dz/2.0;
      const double limitZ  = dz/2.0;

      // Bot outer
      // 0
      out << " 0"
	  << " " << originX << " " << originY << "" << originZ << std::endl;

      // 1
      out << " 1"
	  << " " << limitX << " " << originY << "" << originZ << std::endl;

      // 2
      out << " 2"
	  << " " << limitX << " " << limitY << "" << originZ << std::endl;

      // 3
      out << " 3"
	  << " " << originX << " " << limitY << "" << originZ << std::endl;

      
      // Bot inner
      // 4
      out << " 4"
	  << " " << originInXBot << " " << originInYBot << "" << originZ << std::endl;

      // 5
      out << " 5"
	  << " " << limitInXBot << " " << originInYBot << "" << originZ << std::endl;

      // 6
      out << " 6"
	  << " " << limitInXBot << " " << limitInYBot << "" << originZ << std::endl;

      // 7
      out << " 7"
	  << " " << originInXBot << " " << limitInYBot << "" << originZ << std::endl;

      // Top outer
      // 8
      out << " 8"
	  << " " << originX << " " << originY << "" << limitZ << std::endl;

      // 9
      out << " 9"
	  << " " << limitX << " " << originY << "" << limitZ << std::endl;

      // 10
      out << " 10"
	  << " " << limitX << " " << limitY << "" << limitZ << std::endl;

      // 11
      out << " 11"
	  << " " << originX << " " << limitY << "" << limitZ << std::endl;

      
      // Top inner
      // 12
      out << " 12"
	  << " " << originInXTop << " " << originInYTop << "" << limitZ << std::endl;

      // 13
      out << " 13"
	  << " " << limitInXTop << " " << originInYTop << "" << limitZ << std::endl;

      // 14
      out << " 14"
	  << " " << limitInXTop << " " << limitInYTop << "" << limitZ << std::endl;

      // 15
      out << " 15"
	  << " " << originInXTop << " " << limitInYTop << "" << limitZ << std::endl;


      // ** Construct triangles
      // take into account the triangle orientation to avoid erroneous face culling
      out << "# FACES(triangles)" << std::endl;

      // * Bot faces

      //  - Bot (-y) vertex 0,1,5,4

      // Triangle 0-5-1
      out << " 0 5 1" << std::endl;
      
      // Triangle 0-4-5
      out << " 0 4 5" << std::endl;

      //  - Right (+x) vertex 1,2,6,5

      // Triangle 1-6-2
      out << " 1 6 2" << std::endl;
      
      // Triangle 1-5-6
      out << " 1 5 6" << std::endl;

      //  - Top (+y) vertex 2,3,7,6

      // Triangle 2-7-3
      out << " 2 7 3" << std::endl;
      
      // Triangle 2-6-7
      out << " 2 6 7" << std::endl;

      //  - Left (-x) vertex 3,0,4,7

      // Triangle 3-7-4
      out << " 3 7 4" << std::endl;
      
      // Triangle 3-4-0
      out << " 3 4 0" << std::endl;


      // * Top faces

      //  - Bot (-y) vertex 8,9,13,12

      // Triangle 8-9-13
      out << " 8 9 13" << std::endl;
      
      // Triangle 8-13-12
      out << " 8 13 12" << std::endl;

      //  - Right (+x) vertex 9,10,14,13

      // Triangle 9-10-14
      out << " 9 10 14" << std::endl;
      
      // Triangle 9-14-13
      out << " 9 14 13" << std::endl;

      //  - Top (+y) vertex 10,11,15,14

      // Triangle 10-11-15
      out << " 10 11 15" << std::endl;
      
      // Triangle 10-15-14
      out << " 10 15 14" << std::endl;

      //  - Left (-x) vertex 11,8,12,15

      // Triangle 11-8-12
      out << " 11 8 12" << std::endl;
      
      // Triangle 11-12-15
      out << " 11 12 15" << std::endl;


      // * Outer lateral faces

      //  - Bot (-y) vertex 0,1,9,8

      // Triangle 0-1-9
      out << " 0 1 9" << std::endl;
      
      // Triangle 1-9-8
      out << " 0 9 8" << std::endl;

      //  - Right (+x) vertex 1,2,10,9
      
      // Triangle 1-2-10
      out << " 1 2 10" << std::endl;
      
      // Triangle 1-10-9
      out << " 1 10 9" << std::endl;

      //  - Top (+y) vertex 2,3,11,10
      
      // Triangle 2-3-11
      out << " 2 3 11" << std::endl;
      
      // Triangle 2-11-10
      out << " 2 11 10" << std::endl;

      //  - Left (-x) vertex 3,0,8,11
      
      // Triangle 3-0-8
      out << " 3 0 8" << std::endl;
      
      // Triangle 3-8-11
      out << " 3 8 11" << std::endl;

      
      // * Inner lateral faces

      //  - Left (-x) vertex 4,7,15,12

      // Triangle 4-7-15
      out << " 4 7 15" << std::endl;
      
      // Triangle 4-15-12
      out << " 4 15 12" << std::endl;

      //  - Top (+y) vertex 7,6,14,15

      // Triangle 7-6-14
      out << " 7 6 14" << std::endl;
      
      // Triangle 7-14-15
      out << " 7 14 15" << std::endl;

      //  - Right (+x) vertex 6,5,13,14

      // Triangle 6-5-13
      out << " 6 5 13" << std::endl;
      
      // Triangle 6-13-14
      out << " 6 13 14" << std::endl;

      //  - Bot (-y) vertex 5,4,12,13

      // Triangle 5-4-12
      out << " 5 4 12" << std::endl;
      
      // Triangle 5-12-13
      out << " 5 12 13" << std::endl;
      
    }

  };
};
