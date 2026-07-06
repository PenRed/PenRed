
//
//
//    Copyright (C) 2019-2021 Universitat de València - UV
//    Copyright (C) 2019-2021 Universitat Politècnica de València - UPV
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
//        vicent.gimenez.alventosa@gmail.com
//        vicente.gimenez@uv.es
//    
//


#ifndef __PEN_GEOMETRY__
#define __PEN_GEOMETRY__


#include "geometry_classes.hh"
#include <algorithm>
#include <functional>
#include <utility>

struct pen_surfDS{
  double S;
  unsigned IS;

  inline void set(const double inS,
		  const int inIS){
    S = inS;
    IS = inIS;
  }

  inline bool operator<(const pen_surfDS& r){
    return (S < r.S);
  }
  inline bool operator>(const pen_surfDS& r){
    return (S > r.S);
  }
};

inline bool operator<(const pen_surfDS& l, const pen_surfDS& r){
  return (l.S < r.S);
}
inline bool operator>(const pen_surfDS& l, const pen_surfDS& r){
  return (l.S > r.S);
}  

struct pen_quadSurface{

  double AXX;
  double AXY;
  double AXZ;
  double AYY;
  double AYZ;
  double AZZ;
  double AX;
  double AY;
  double AZ;
  double A0;

  unsigned int KPLANE; //Control if this surface is a plane
  
};

struct pen_bodySurf{
  unsigned int  KSURF;
  unsigned int  KFLAG;

  inline void set(const unsigned ksurf, const unsigned kflag){
    KSURF = ksurf;
    KFLAG = kflag;
  }

  inline void get(unsigned& ksurf, unsigned& kflag) const {
    ksurf = KSURF;
    kflag = KFLAG;
  }
  
};

struct pen_quadBody : public pen_baseBody{

  static const unsigned int NXG = 250;   //Maximum number of neighbors ??
  
  char BALIAS[5];

  unsigned int KBODY[NXG];
  unsigned int KBOMO;

  unsigned int  KMOTH;
  unsigned int  KDGHT[NXG];

  pen_bodySurf surfs[NXG];
  
  inline void setSurf(const unsigned index,
		      const unsigned ksurf,
		      const unsigned kflag){
    surfs[index].set(ksurf,kflag);
  }
  inline void getSurf(const unsigned index,
		      unsigned& ksurf,
		      unsigned& kflag) const{
    surfs[index].get(ksurf,kflag);
  }
  
  inline void setKSURF(const unsigned index, const unsigned ksurf){
    surfs[index].KSURF = ksurf;
  }
  inline void setKFLAG(const unsigned index, const unsigned kflag){
    surfs[index].KFLAG = kflag;
  }
  
  inline unsigned getKSURF(const unsigned index) const{
    return surfs[index].KSURF;
  }
  inline unsigned getKFLAG(const unsigned index) const{
    return surfs[index].KFLAG;
  }  
};

class pen_quadricGeo : public abc_geometry<pen_quadBody>{
  DECLARE_GEOMETRY(pen_quadricGeo)
protected:
  
  const double FUZZL = 1.0E-12;  
  const double mFUZZL = -FUZZL;  
  const double FUZZT = -0.25E0*FUZZL;
  static const unsigned int NS  = 10000; //Maximum number of surfaces
  static const unsigned int NS2M = 2*NS;
  //Surfaces stack
  unsigned int NSURF;
  pen_quadSurface surfaces[NS];
  
  bool LVERB;
  unsigned NMATG;

  virtual penred::errors::Error GEOMIN(FILE* IRD,
				       FILE* IWR,
				       const unsigned verbose);

  void move(const double DS, pen_particleState& state) const;
  bool inBody(const pen_quadBody* pbody, const unsigned KSP[NS]) const;
  bool inBody(const unsigned KB, const unsigned KSP[NS]) const;
  bool goInner(const pen_quadBody*& pbody,
	       pen_particleState& state,
	       unsigned KSP[NS],
	       pen_surfDS surfDS[NS2M],
	       int& NSC,
	       unsigned& KSLAST) const;
  
  void STEPSI(const pen_quadBody* pbody, const pen_particleState& state, unsigned KSP[NS], pen_surfDS surfDS[NS2M], int &NSC_IO, unsigned& KSLAST_IO) const;
  void STEPLB(const pen_quadBody* pbody, pen_particleState& state, const unsigned KSP[NS], int &IERR) const;
  
public:

  enum errors{
    SUCCESS = 0,
    MISSING_PARAMETER,
    BAD_VALUE,
    UNABLE_TO_OPEN_FILE,
    NO_FILE_PROVIDED,
    GEOMIN_ERROR,
    NB_LIMIT_REACHED,
    NS_LIMIT_REACHED,
    NXG_LIMIT_REACHED,
    UNRESOLVED_BODY,
    INCONSISTENT_LABEL,
    UNKNOWN_PARTICLE,
    UNKNOWN_ERROR,
  };

  static constexpr const char* errorMessage(const int val) noexcept {
    switch(val){
    case SUCCESS: return "Success";
    case MISSING_PARAMETER: return "Missing parameter";
    case BAD_VALUE: return "Invalid value";
    case UNABLE_TO_OPEN_FILE: return "Unable to open file";
    case NO_FILE_PROVIDED: return "No file provided";
    case GEOMIN_ERROR: return "GEOMIN rutine error";
    case NB_LIMIT_REACHED: return "Maximum number of bodies reached";
    case NS_LIMIT_REACHED: return "Maximum number of surfaces reached";
    case NXG_LIMIT_REACHED: return "Maximum number of surfaces per body reached";
    case UNKNOWN_PARTICLE: return "Unknown particle";
    case UNRESOLVED_BODY: return "Unable to resolve body";
    default: return "Unknown error";
    }
  }

  static const unsigned int NXG = pen_quadBody::NXG;
  
  pen_quadricGeo();

  inline const pen_quadSurface& getSurface(unsigned int ks) const {
    if(ks >= NSURF){
      char error[300];
      sprintf(error,"getSurface: Ordered surface (%u) out of range (%d)",ks,NSURF);
      throw std::out_of_range(error);	      
    }
    return surfaces[ks];
  }
  inline unsigned int getSurfaces() const {return NSURF;}
  unsigned getIBody(const char* elementName) const final override;
  std::string getBodyName(const unsigned ibody) const final override;
  
  
  penred::errors::Error specificConfigure(const pen_parserSection& config,
					  const unsigned verbose) final override;
  void locateLocal(pen_particleState& state) const final override;
  void stepLocal(pen_particleState& state,
                 double DS,
                 double &DSEF,
                 double &DSTOT,
                 int &NCROSS) const final override;
  //PENGEOM_mod (except DSTOT and KSLAST. This varaibles will be passed as STEP arguments)
  //QSURF
  
};

inline void pen_quadricGeo::move(const double DS, pen_particleState& state) const{
  state.X += DS*state.U;
  state.Y += DS*state.V;
  state.Z += DS*state.W;
}

inline bool pen_quadricGeo::inBody(const pen_quadBody* pbody,
				   const unsigned KSP[NS]) const{

  // Checks if the particle with the KSP side pointers is in the
  // the specified body i.e. if the particle has been
  // crossed a limiting surface (KFLAG < 3) of KB.
  
  // pbody  -> Body pointer
  // KSP    -> Body surface side pointers respect the particle direction

  const pen_bodySurf* psurfLast = pbody->surfs+pbody->getKSURF(NXG-1);
  for(const pen_bodySurf* psurf = pbody->surfs; psurf != psurfLast; ++psurf)
    {
      if(psurf->KFLAG < 3 && KSP[psurf->KSURF-1] != psurf->KFLAG){
	//Limiting surface has been crossed, is not in
	//that body
	return false;
      }
    }

  return true;
}

inline bool pen_quadricGeo::inBody(const unsigned KB, const unsigned KSP[NS]) const{
  return inBody(bodies + (KB-1),KSP);
}

inline bool pen_quadricGeo::goInner(const pen_quadBody*& pbody,
				    pen_particleState& state,
				    unsigned KSP[NS],
				    pen_surfDS surfDS[NS2M],
				    int& NSC,
				    unsigned& KSLAST) const{

  //Go inner in the geometry tree until a material body has been found.
  //As input, the initial body (KB) and the NSC surface side pointers (KSP)
  //must be specified. Use STEPSI to obtain the surface
  //side pointers.

  //As output, the final body containing the particle will be set to
  //KB. Also, the particle  IBODY and MAT  will be updated. KSP,
  //surfDS, NSC and KLAST will be updated by the last call to STEPSI.

  //The function returns "true" if the particle is outside the input
  //body or module
  
  //Try to down one step in geometry tree from input bodt 
  int IERR;
  STEPLB(pbody,state,KSP,IERR);
  //Update body index
  pbody = bodies + (state.IBODY-1);
  
  //Iterate until the particle reaches a material body
  while(IERR != 0){
    if(IERR == -1)
      {
	//  ****  The particle enters a submodule.
	STEPSI(pbody,state,KSP,surfDS,NSC,KSLAST);
	STEPLB(pbody,state,KSP,IERR);
	//Update body index
	pbody = bodies + (state.IBODY-1);
      }
    else if(IERR == 1){
      //  ****  The particle is not in the body or module.
      if(state.IBODY <= NBODYS)
	{
	  //Is in another body, update side pointers and go inner
	  STEPSI(pbody,state,KSP,surfDS,NSC,KSLAST);
	  STEPLB(pbody,state,KSP,IERR);
	  //Update body index
	  pbody = bodies + (state.IBODY-1);
	}
      else
	{
	  //  ****  The particle leaves the enclosure.
	  return true;
	}
    }
  }
  return false;
  
  
}

inline void FSURF(const pen_quadSurface& surface,
		  const pen_particleState& state,
		  double &A,
		  double &B,
		  double &C)
{
  //     Calculates the parameters of the master function of the surface KS
  //  and the ray (X,Y,Z)+S*(U,V,W).


  if(surface.KPLANE == 0)
    {
      A = state.U*(surface.AXX*state.U+
		      surface.AXY*state.V+
		      surface.AXZ*state.W)+
	state.V*(surface.AYY*state.V+
		    surface.AYZ*state.W)+
	state.W*surface.AZZ*state.W;
      
      double XXX = surface.AXX*state.X+surface.AXY*state.Y+
	surface.AXZ*state.Z+surface.AX;
      
      double YYY = surface.AYY*state.Y+surface.AYZ*state.Z+surface.AY;
      double ZZZ = surface.AZZ*state.Z+surface.AZ;

      B = state.U*(surface.AXX*state.X+XXX)+
	state.V*(surface.AXY*state.X+surface.AYY*state.Y+YYY)+
	state.W*(surface.AXZ*state.X+surface.AYZ*state.Y+surface.AZZ*state.Z+
		    ZZZ);

      C = state.X*XXX+state.Y*YYY+state.Z*ZZZ+surface.A0;
    }
  else
    {
      A = 0.0;
      B = state.U*surface.AX+state.V*surface.AY+state.W*surface.AZ;
      C = state.X*surface.AX+state.Y*surface.AY+state.Z*surface.AZ+surface.A0;
    }
}

void ROTSHF(double OMEGA,
	    double THETA,
	    double PHI,
	    double DX,
	    double DY,
	    double DZ,
	    double &AXX,
	    double &AXY,
	    double &AXZ,
	    double &AYY,
	    double &AYZ,
	    double &AZZ,
	    double &AX,
	    double &AY,
	    double &AZ,
	    double &A0);

template<>
constexpr const char* penred::errors::errorMessage<pen_quadricGeo>(const int val) noexcept {
  return pen_quadricGeo::errorMessage(val);
}

#endif
