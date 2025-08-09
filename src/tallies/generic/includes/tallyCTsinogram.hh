
//
//
//    Copyright (C) 2020-2021 Universitat de València - UV
//    Copyright (C) 2020-2021 Universitat Politècnica de València - UPV
//    Copyright (C) 2025 Vicent Giménez Alventosa
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
  DECLARE_TALLY(pen_CTsinogram,pen_particleState,CT_SINOGRAM,
		std::pair<double, penred::tally::Dim<2>>)
  
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
				      USE_SAMPLEDPART |
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
  {
    setResultsGenerator<0>
      ([this](const unsigned long long nhists) -> penred::measurements::results<double, 2>{
	
	double invn = 1.0/static_cast<double>(nhists);
	    
	//Create results
	penred::measurements::results<double, 2> results;
	results.initFromLists({static_cast<unsigned>(npixs), nphi},
			      {penred::measurements::limitsType(0.0, phi*180.0/M_PI),
			       penred::measurements::limitsType(phi0*180.0/M_PI, phif*180.0/M_PI)});

	results.description =
	  "PenRed: CT Sinogram Report\n\n"
	  "  This tally estimates the product of the linear attenuation coefficient (\u03BC)\n"
	  "  and the path length (x) through the object (\u03BC\u00B7x) for each detector bin\n"
	  "  and projection.\n"
	  "  The value is derived from the Beer-Lambert law: \u03BC\u00B7x = ln(N\u2080/N),\n"
	  "  where N\u2080 is the fluence at the source and N is the detected fluence.\n"
	  "  The detector geometry is simplified as an arc, with angular positions for each\n"
	  "  bin reported in degrees.\n\n";
	results.description += " Histories simulated: " + std::to_string(nhists) + "\n\n";

	results.setDimHeader(0, "Detection angle (deg)");
	results.setDimHeader(1, "Projection Angle (deg)");
	results.setValueHeader("  \u03BC\u00B7x  ");

	for(size_t j = 0; j < sinoDim; ++j){
	  double q, q2,sigma;
	  double qNorm, q2Norm,sigmaNorm;
	  double sigmaMux;
            
	  q = sino[j]*invn;
	  q2 = sino2[j]*invn;
	  sigma = (q2-(q*q))*invn;
	  if(sigma > 0.0){ sigma = sqrt(sigma);}
	  else{sigma = 0.0;}
            
	  qNorm = sinoNorm[j]*invn;
	  q2Norm = sino2Norm[j]*invn;
	  sigmaNorm = (q2Norm-(qNorm*qNorm))*invn;
	  if(sigmaNorm > 0.0){ sigmaNorm = sqrt(sigmaNorm);}
	  else{sigmaNorm = 0.0;}

	  //Estimate mu*x as ln(qNorm/q)
	  //
	  //  N=N0*e^(-mu*x) -> ln(N/N0) = -mu*x -> ln(N0/N) = mu*x
	  //
	  double mux = 0.0E0;
	  if(q > 0.0E0){
	    mux = qNorm/q;
	    if(mux > 0.0){
	      mux = log(mux);
	    } else{
	      mux = 0.0;
	    }   
	    sigmaMux = sqrt(pow(sigmaNorm/qNorm,2) + pow(sigma/q,2));
	  }
	  else{
	    sigmaMux = 0.0;
	  }

	  results.data[j] = mux;
	  results.sigma[j] = sigmaMux;
	}

	return results;
      });    
  }
    
        
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
  void tally_sampledPart(const unsigned long long nhist,
			 const unsigned long long dhist,
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
