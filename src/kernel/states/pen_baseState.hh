
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



#ifndef __PENRED_BASIC_PARTICLE_STATE_
#define __PENRED_BASIC_PARTICLE_STATE_

#include <cmath>
#include <string>

//-------------------
// Particle basic state
//-------------------

inline const char* baseStateHeader(){
  return "     E(eV)        X(cm)        Y(cm)        Z(cm)         U            V            W           WGHT     IBODY  MAT      PAGE       ILBs";
}

inline const char* baseStateHeaderNoGeo(){
  return "     E(eV)        X(cm)        Y(cm)        Z(cm)         U            V            W           WGHT          PAGE       ILBs";
}

struct pen_particleState{

protected:

public:

    static const size_t ILBsize = 5*sizeof(int);

  //  ****  Particle TRACK variables (to be initialised before calling
  //        subroutine START).  
  //  ----  Energy, position, direction, and weight.  
  double E, X, Y, Z, U, V, W, WGHT;
  //  ----  Particle type, current body, and material.
  unsigned IBODY, MAT;
  //  ----  Particle history flags.
  int ILB[5];
  //  ****  The particle age (time elapsed since the start of the shower)
  //        is recorded when LAGE=.TRUE.
  bool LAGE;
  double PAGE;

pen_particleState() : E(0.0), X(0.0), Y(0.0), Z(0.0),
    U(0.0), V(0.0), W(1.0), WGHT(1.0),
    IBODY(0), MAT(0), LAGE(false), PAGE(0.0)
  {
    ILB[0] = 0;
    ILB[1] = 0;
    ILB[2] = 0;
    ILB[3] = 0;
    ILB[4] = 0;    
  }

pen_particleState(const pen_particleState& c){
  *this = c;
}

  inline double posMod() const {return sqrt(X*X+Y*Y+Z*Z);}
  inline double posMod2() const {return X*X+Y*Y+Z*Z;}
  
  void baseReset(){
    E = 0.0;
    X = 0.0;
    Y = 0.0;
    Z = 0.0;
    U = 0.0;
    V = 0.0;
    W = 1.0;
    WGHT = 1.0;
    IBODY = 0;
    MAT = 0;
    PAGE = 0.0;

    memset(ILB,0,ILBsize);
  }
  virtual void reset(){
    E = 0.0;
    X = 0.0;
    Y = 0.0;
    Z = 0.0;
    U = 0.0;
    V = 0.0;
    W = 1.0;
    WGHT = 1.0;
    IBODY = 0;
    MAT = 0;
    PAGE = 0.0;

    memset(ILB,0,ILBsize);
  }

  //Ensure that input variable is a derived class from "pen_particleState"
  /*
    
    const typename std::enable_if<
    std::is_base_of<pen_particleState, derivedState>::value,
    derivedState&
    >::type
  */
  
  inline void copyBase(const pen_particleState& cstate){
    *this = cstate;
  }
    
  inline std::string stringifyBase() const {

    char text[300];
    sprintf(text,"%12.4E %12.4E %12.4E %12.4E "
	    "%12.4E %12.4E %12.4E "
	    "%12.4E  %4u  %4u %12.4E "
	    "  %d %d %d %d %d",
	    E,X,Y,Z,
	    U,V,W,
	    WGHT,IBODY,MAT,PAGE,
	    ILB[0],ILB[1],ILB[2],ILB[3],ILB[4]);

    return std::string(text);
  }

  inline std::string stringifyBaseNoGeo() const {

    char text[300];
    sprintf(text,"%12.4E %12.4E %12.4E %12.4E "
	    "%12.4E %12.4E %12.4E "
	    "%12.4E  %12.4E "
	    "  %d %d %d %d %d",
	    E,X,Y,Z,
	    U,V,W,
	    WGHT,PAGE,
	    ILB[0],ILB[1],ILB[2],ILB[3],ILB[4]);

    return std::string(text);
  }
  
  virtual std::string stringify() const {
    return stringifyBase();
  }
  
};

//Only for debug
inline bool checkNANState(const pen_particleState& state, const unsigned verbose = 0){

  if(std::isnan(state.E) || std::isnan(state.X) || std::isnan(state.Y) ||
     std::isnan(state.Z) || std::isnan(state.U) || std::isnan(state.V) ||
     std::isnan(state.W) || std::isnan(state.WGHT) ||
     std::isnan(state.PAGE)){

    if(verbose > 0)
      printf("%s\n",state.stringifyBase().c_str());
    return true;
  }
  return false;
}


//Create a generic copy function

// Case 1: Both state types are equal
template<class Tstate1, class Tstate2>

// Ensure Tstate1 == Tstate2 using return type (void)
inline typename std::enable_if<std::is_same<Tstate1,Tstate2>::value,void>::type  

stateCopy(Tstate1& state1, const Tstate2& state2){
  state1 = state2;
}

// Case 2: States types are different but state1 is convertible to state2
template<class Tstate1, class Tstate2>

// Ensure Tstate1 != Tstate2 using return type (void)
inline typename std::enable_if<!(std::is_same<Tstate1,Tstate2>::value),void>::type  

stateCopy(Tstate1& state1, const Tstate2& state2,
	  // Ensure Tstate1 is convertible to Tstate2 using a extra dummy argument
	  const typename std::enable_if<std::is_convertible<Tstate1&,Tstate2&>::value,void>::type* /*dummy*/ = nullptr)
{
  //Get a reference to Tstate2 object from state1
  Tstate2& conv = static_cast<Tstate2&>(state1);

  //Fill Tstate2 object part of state1 with state2 object
  conv = state2;
}

// Case 3: States types are different but state2 is convertible to state1
template<class Tstate1, class Tstate2>

// Ensure Tstate1 != Tstate2 using return type (void)
inline typename std::enable_if<!(std::is_same<Tstate1,Tstate2>::value),void>::type 

stateCopy(Tstate1& state1, const Tstate2& state2,
	  // Ensure Tstate2 is convertible to Tstate1 using a extra dummy argument
	  const typename std::enable_if<std::is_convertible<Tstate2&,Tstate1&>::value,void>::type* /*dummy*/ = nullptr)
{
  //Get a reference to Tstate1 object from state2
  const Tstate1& conv = static_cast<const Tstate1&>(state2);

  //Fill Tstate2 object part of state1 with state2 object
  state1 = conv;
}

// Case 4: State types are different and state1 (state2) isn't convertible to state2 (state1)
template<class Tstate1, class Tstate2>

// Ensure Tstate1 != Tstate2 using return type (void)
inline typename std::enable_if<!(std::is_same<Tstate1,Tstate2>::value),void>::type  

stateCopy(Tstate1& state1, const Tstate2& state2,
	  // Ensure Tstate1 is not convertible to Tstate2 using a extra dummy argument
	  const typename std::enable_if<!(std::is_convertible<Tstate1&,Tstate2&>::value),void>::type* /*dummy12*/ = nullptr,
	  // Ensure Tstate2 is not convertible to Tstate1 using a extra dummy argument
	  const typename std::enable_if<!(std::is_convertible<Tstate2&,Tstate1&>::value),void>::type* /*dummy21*/ = nullptr)
{

  //Copy only base state
  state1.copyBase(state2);
}

#endif
