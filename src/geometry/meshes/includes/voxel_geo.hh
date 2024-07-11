
//
//
//    Copyright (C) 2019-2023 Universitat de València - UV
//    Copyright (C) 2019-2023 Universitat Politècnica de València - UPV
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


#ifndef __PENRED_VOXELIZED_GEOMETRY__
#define __PENRED_VOXELIZED_GEOMETRY__

#include <cmath>
#include <stdlib.h>
#include "pen_dump.hh"

#ifdef _PEN_USE_DICOM_
#include "pen_dicom.hh"
#endif

struct pen_voxel : pen_baseMesh{

  double densityFact; //density(voxel)/density(material)

  pen_voxel() : densityFact(1.0)
  {}
};

class pen_voxelGeo : public abc_mesh<pen_voxel>{
  DECLARE_GEOMETRY(pen_voxelGeo)
  
  protected:
  //Number of voxels
  unsigned nx, ny, nz;
  unsigned long nxy;
  //Voxels size
  double dx, dy, dz;
  double idx, idy, idz;
  //Geometry limits (low limits are 0,0,0 for each axy)
  double Mdx, Mdy, Mdz;
  
  //Enclosure variables
  double enclosureMargin;
  double enclosureXlimit, enclosureYlimit, enclosureZlimit;
  double enclosureXlimit0, enclosureYlimit0, enclosureZlimit0;

  void move(const double ds, pen_particleState& state) const;
  bool crossVox(const double ds,
		const unsigned imat,
		const unsigned nvoxAxis,
		const int dvox1D,
		const long int dvox3D,
		double& remaining,
		unsigned long& ivox,
		long int& iaxis,
		double& DSEF,
		double& DSTOT,
		int& NCROSS,
		pen_particleState& state) const;
  
public:

  // Define "zero" value to handle particles situated
  // on voxels and mesh limits 
  static const double pZero;
  static const double nZero;
  static const double inf;
  
  pen_voxelGeo();
  
  inline unsigned xVox() const {return nx;}
  inline unsigned yVox() const {return ny;}
  inline unsigned zVox() const {return nz;}
  inline unsigned nVox() const {return nElements;}
  
  inline double xSize() const {return dx;}
  inline double ySize() const {return dy;}
  inline double zSize() const {return dz;}

  
  
  virtual int configure(const pen_parserSection& config, const unsigned verbose) override;

  void locate(pen_particleState& state) const final override;
  void locateInMesh(pen_particleState& state) const;

  void step(pen_particleState& state,
	    double DS,
	    double &DSEF,
	    double &DSTOT,
	    int &NCROSS) const final override;
  void stepInMesh(pen_particleState& state,
	    double DS,
	    double &DSEF,
	    double &DSTOT,
	    int &NCROSS) const;
        
  bool enterEnclosure(const double x, const double y, const double z,
                      const double u, const double v, const double w,
                      double &ds) const;

  void exitEnclosure(const double x, const double y, const double z,
		     const double u, const double v, const double w,
		     double &ds) const;
        
  bool crossMesh(const double x,
                 const double y,
                 const double z,
                 const double u,
                 const double v,
                 const double w,
                 double& ds) const;  
  
  int setVoxels(const unsigned nvox[3],
		const double size[3],
		const unsigned* mats,
		const double* dens,
		const unsigned verbose = 0);

  int loadData(const unsigned char* data,
	       size_t& pos,
	       const unsigned verbose = 0);

  int dump(unsigned char*& data,
	   size_t& pos,
	   const unsigned verbose) const;

  int dump2File(const char* filename);
  
  int loadFile(const char* filename,
	       const unsigned verbose = 0);

  int loadASCII(const char* filename,
		const unsigned verbose = 0);  

  unsigned getIBody(const char* bname) const override;

  std::string getBodyName(const unsigned ibody) const override;
  
  virtual int printImage(const char* filename) const;
};

inline void pen_voxelGeo::move(const double ds, pen_particleState& state) const{
  state.X += ds*state.U;
  state.Y += ds*state.V;
  state.Z += ds*state.W;			    
}

inline bool pen_voxelGeo::crossVox(const double ds,
				   const unsigned imat,
				   const unsigned nvoxAxis,
				   const int dvox1D,
				   const long int dvox3D,
				   double& remaining,
				   unsigned long& ivox,
				   long int& iaxis,
				   double& DSEF,
				   double& DSTOT,
				   int& NCROSS,
				   pen_particleState& state) const {

  //
  // ds -> distance to next wall
  // imat -> initial material index
  // nvoxAxis -> number of voxels in the axis that will be crossed
  // dvox1D -> Index increment in crossed axis (-1 or +1 for backward or forward)
  // dvox3D -> Global index increment (+-1 for X, +- nx for Y and +- nx*ny for Z)
  // remaining -> Remaining distance to be traveled by the particle
  // ivox -> Actual global (3D) voxel index
  // iaxis -> Actual voxel index in crossed axis (1D) (ix,iy,iz)
  // DSEF -> Total traveled distance in original material
  // DSTOT -> Total traveled distance
  // NCROSS -> Number of crossed surfaces
  // state -> Actual particle state
  //
  
  if(ds > 0.0E0){
    // Calculate distance consumed at voxel density
    double consumed = ds*mesh[ivox].densityFact;
    if(consumed > remaining){  // Maximum traveled distance reached
      DSEF += remaining/mesh[ivox].densityFact;
      DSTOT = DSEF;
      // Update IBODY and material indexes
      state.IBODY = imat;
      state.MAT = imat;

      // Update position
      move(DSEF,state);

      return true;
    }
    DSEF += ds;    // Update traveled distance
    // Update remaining traveling distance into material density
    remaining -= consumed;
	
  }

  iaxis += dvox1D;  // Update axis voxel index
  if(iaxis < 0 || iaxis >=(long int)nvoxAxis){ // Particle scapes from geometry	
      
    // Crossing the enclosure counts as an interface, regardless the
    // material of the enclosure
    NCROSS = 1;
    
    state.IBODY = 0; //Go to enclosure
    state.MAT = enclosureMat;
    DSTOT = DSEF; //Particle scape to the enclosure

    // Update position
    move(DSEF,state);
	
    return true;
  }

  // Update global voxel index
  ivox += dvox3D;
  unsigned vmat = mesh[ivox].MATER;
  if(imat != vmat){ //Material changed, finish step	
    NCROSS = 1;
    state.IBODY = vmat;
    state.MAT = vmat;
    DSTOT = DSEF;
    
    // Update position
    move(DSEF,state);

    return true;
  }

  return false;
}

bool moveIn(double pos,
	    const double dir,
	    double& ds2in,
	    double& ds2out,
	    const double max);
  
inline bool pen_voxelGeo::crossMesh(const double x, 
                                    const double y, 
                                    const double z,
                                    const double u,
                                    const double v, 
                                    const double w,
                                    double& ds) const{
        
    //Check if the particle reaches the mesh in
    //each axis
    double dsxIn, dsyIn, dszIn;
    double dsxOut, dsyOut, dszOut;
    if(!moveIn(x, u, dsxIn, dsxOut, Mdx) ||
       !moveIn(y, v, dsyIn, dsyOut, Mdy) ||
       !moveIn(z, w, dszIn, dszOut, Mdz)){
    
      return false;
    }

    //Move the particle the required distance to reach
    //geometry mesh on all axis.
    double ds2allIn = std::max(std::max(dsxIn,dsyIn),dszIn);
    //Check if the particle go out of the mesh when this distance is traveled
    if(ds2allIn > 0.0){
      double ds2someOut = std::min(std::min(dsxOut,dsyOut),dszOut);
      if(ds2allIn >= ds2someOut){
        return false;
      }
    }
    
    ds = ds2allIn;
    
    return true;
}

inline bool pen_voxelGeo::enterEnclosure(const double x, 
					 const double y, 
					 const double z,
					 const double u,
					 const double v, 
					 const double w,
					 double& ds) const{
        
  //Check if the particle reaches the mesh in
  //each axis
  double dsxIn, dsyIn, dszIn;
  double dsxOut, dsyOut, dszOut;
  if(!moveIn(x+enclosureMargin, u, dsxIn, dsxOut, enclosureXlimit0) ||
     !moveIn(y+enclosureMargin, v, dsyIn, dsyOut, enclosureYlimit0) ||
     !moveIn(z+enclosureMargin, w, dszIn, dszOut, enclosureZlimit0)){
    
    return false;
  }

  //Move the particle the required distance to reach
  //geometry mesh on all axis.
  double ds2allIn = std::max(std::max(dsxIn,dsyIn),dszIn);
  //Check if the particle go out of the mesh when this distance is traveled
  if(ds2allIn > 0.0){
    double ds2someOut = std::min(std::min(dsxOut,dsyOut),dszOut);
    if(ds2allIn >= ds2someOut){
      return false;
    }
  }
    
  ds = ds2allIn;
    
  return true;
}

inline void pen_voxelGeo::exitEnclosure(const double x, 
					const double y, 
					const double z,
					const double u,
					const double v, 
					const double w,
					double& ds) const{
        
  //Check if the particle reaches the mesh in
  //each axis
  double dsxIn, dsyIn, dszIn;
  double dsxOut, dsyOut, dszOut;
  moveIn(x+enclosureMargin, u, dsxIn, dsxOut, enclosureXlimit0);
  moveIn(y+enclosureMargin, v, dsyIn, dsyOut, enclosureYlimit0);
  moveIn(z+enclosureMargin, w, dszIn, dszOut, enclosureZlimit0);

  //Move the particle the required distance to reach
  //geometry mesh on some axis.
  ds = std::min(std::min(dsxOut,dsyOut),dszOut);
}

#endif
