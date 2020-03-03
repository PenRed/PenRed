
//
//
//    Copyright (C) 2019 Universitat de València - UV
//    Copyright (C) 2019 Universitat Politècnica de València - UPV
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
  unsigned nx, ny, nz, nxy;
  //Voxels size
  double dx, dy, dz;
  double idx, idy, idz;
  //Geometry limits (low limits are 0,0,0 for each axy)
  double Mdx, Mdy, Mdz;

  void move(const double ds, pen_particleState& state) const;
  bool crossVox(const double ds,
		const unsigned imat,
		const unsigned nvoxAxis,
		const int dvox1D,
		const int dvox3D,
		double& remaining,
		size_t& ivox,
		int& iaxis,
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

  
  
  virtual int configure(const pen_parserSection& config, const unsigned verbose);

  void locate(pen_particleState& state) const final override;

  void step(pen_particleState& state,
	    double DS,
	    double &DSEF,
	    double &DSTOT,
	    int &NCROSS) const final override;
  
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

  unsigned getIBody(const char* bname) const;
  
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
				   const int dvox3D,
				   double& remaining,
				   size_t& ivox,
				   int& iaxis,
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
      state.IBODY = imat-1;
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
  if(iaxis < 0 || iaxis >= (int)nvoxAxis){ // Particle scapes from geometry	
    NCROSS = 1;
    state.IBODY = constants::MAXMAT;
    state.MAT = 0;
    DSTOT = inf; //Particle scape, move it to "inf"

    // Update position
    move(inf,state);
	
    return true;
  }

  // Update global voxel index
  ivox += dvox3D;
  unsigned vmat = mesh[ivox].MATER;
  if(imat != vmat){ //Material changed, finish step	
    NCROSS = 1;
    state.IBODY = vmat-1;
    state.MAT = vmat;
    DSTOT = DSEF;
    
    // Update position
    move(DSEF,state);

    return true;
  }

  return false;
}


bool moveIn(const double pos,
	    const double dir,
	    double& ds,
	    const double max);
#endif
