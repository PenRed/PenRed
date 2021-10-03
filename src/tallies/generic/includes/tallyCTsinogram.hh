
//
//
//    Copyright (C) 2020-2021 Universitat de València - UV
//    Copyright (C) 2020-2021 Universitat Politècnica de València - UPV
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
//        sanolgi@upvnet.upv.es (Sandra Oliver Gil)
//        vicent.gimenez.alventosa@gmail.com  (Vicent Giménez Alventosa)
//        vicente.gimenez@uv.es (Vicente Giménez Gómez)
//    
//

#ifndef __PEN_CT_SINOGRAM_TALLY__
#define __PEN_CT_SINOGRAM_TALLY__

#include "pen_constants.hh" 

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


struct CTsinogram_point
  {
   double X,Y,Z;
  };

class pen_CTsinogram : public pen_genericTally<pen_particleState> {
  DECLARE_TALLY(pen_CTsinogram,pen_particleState)
  
  private:
  //Position of the center of the virtual ring
  double xOrigin, yOrigin, zOrigin;
  //Value of the detector/pixel height and its semi value
  double dz, semiDet;
  //Variables of the radious of the virtual ring, inner radious ri,
  //outer radious ro.
  //Variables of the angular information, initial anf final angle,
  //phi0 phif,angular step dphi, inner arc of the detector si,
  //angular size of the detector phi and number of angular pixels npixs
  double  ri, ri2,// ro,
    phi0, phif, dphi,
    si, phi, phipx;
  //Time information, minimum and maximum time and time window
  double tmin,tmax,dt;
  //Energy interval
  double emin, emax;
  //Number of pixels of the detector
  int npixs;
  //Number of angular projections
  unsigned long nphi;
  //Maximum number of pixels of the detector
  static const int npxmax=32000;
  //Kind of particle
  pen_KPAR part;
  //Sinogram dimension and sinogram memory size
  size_t sinoDim;
  double* sino;
  double* sino2;
  double* sinotemp;
  double* lasthist;

  double* sinoNorm;
  double* sino2Norm;
  double* sinotempNorm;
  double* lasthistNorm;
  
  CTsinogram_point lastPos;
  double lastPage;
  bool inGeo;
  bool scatter = true;
  std::vector<pen_EdepMat> edepMat;
  
  int nmat;
  int iproj;
  bool knocked = false;
  bool moved2geo = false;
public:
    
  pen_CTsinogram() : pen_genericTally(USE_LOCALEDEP |
				      USE_JUMP |
				      USE_BEGINPART |
				      USE_ENDPART |
				      USE_BEGINHIST |
				      USE_STEP |
				      USE_ENDHIST |
				      USE_MOVE2GEO |
				      USE_ENDSIM |
                      USE_KNOCK),
		     sino(nullptr),
		     sino2(nullptr),
		     sinotemp(nullptr),
		     lasthist(nullptr),
		     sinoNorm(nullptr),
		     sino2Norm(nullptr),
		     sinotempNorm(nullptr),
		     lasthistNorm(nullptr)     
  {}
    
        
    bool detIndexes(const pen_particleState& state,
                    CTsinogram_point finalPoint,
                    const double partPage,
                    unsigned long& iphi, 
                    unsigned long& ipix);
    
    void countDetector(const  unsigned long long nhist, 
                       const pen_particleState& state);
    
    void countSource(const  unsigned long long nhist, 
                     const pen_particleState& state);

  void tally_localEdep(const unsigned long long nhist,
		       const pen_KPAR kpar,
		       const pen_particleState& state,
		       const double dE);
  void tally_beginHist(const unsigned long long nhist,
		       const unsigned kdet,
		       const pen_KPAR kpar,
		       const pen_particleState& state);
  void tally_step(const unsigned long long nhist,
		  const pen_KPAR kpar,
		  const pen_particleState& state,
		  const tally_StepData& stepData);
  
  void tally_move2geo(const unsigned long long nhist,
		      const unsigned kdet,
		      const pen_KPAR kpar,
		      const pen_particleState& state,
		      const double dsef,
		      const double dstot);
  void tally_endHist(const unsigned long long nhist);  
  void tally_jump(const unsigned long long /*nhist*/,
				const pen_KPAR kpar,
				const pen_particleState& state,
				const double /*ds*/);
  
  void tally_beginPart(const unsigned long long nhist,
		       const unsigned kdet,
		       const pen_KPAR kpar,
		       const pen_particleState& state);
  
  void tally_endPart(const unsigned long long nhist,
				   const pen_KPAR kpar,
				   const pen_particleState& state);

  void tally_endSim(const unsigned long long /*nhist*/);
  void tally_knock(const unsigned long long /*nhist*/,
		 const pen_KPAR /*kpar*/,
		 const pen_particleState& /*state*/,
		 const int /*icol*/);
  
  int configure(const wrapper_geometry& geometry,
		const abc_material* const materials[constants::MAXMAT],     
		const pen_parserSection& config, const unsigned verbose);
  void flush();
  void saveData(const unsigned long long nhist) const;
  int sumTally(const pen_CTsinogram& tally);

  inline long int getProjection(const pen_particleState& state) const{
      if(state.PAGE < tmin || state.PAGE >= tmax)
        return -1;
      return (state.PAGE-tmin)/dt;
  }
  
  ~pen_CTsinogram(){
    if(sino != nullptr)
        free(sino);
    if(sino2 != nullptr)
        free(sino2);
    if(sinotemp != nullptr)
        free(sinotemp);
    if(lasthist != nullptr)
        free(lasthist);

    if(sinoNorm != nullptr)
        free(sinoNorm);
    if(sino2Norm != nullptr)
        free(sino2Norm);
    if(sinotempNorm != nullptr)
        free(sinotempNorm);
    if(lasthistNorm != nullptr)
        free(lasthistNorm);
    
  }
};


#endif
