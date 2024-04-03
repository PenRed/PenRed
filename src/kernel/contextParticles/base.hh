
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

 
#ifndef __PENRED_CONTEXT_PARTICLES_BASE_INCLUDE__
#define __PENRED_CONTEXT_PARTICLES_BASE_INCLUDE__

#include "../particles/includes/pen_particles.hh"

namespace penred{

  namespace context{

    //  ** Context particle types related functions **
    //
    // This family of functions is intended to use context particles
    // automatically regardless the context type, allowing templated
    // functions to work for all defined context and particles.

    // ** Particle types definition
    
    // ** Define a structure to get the tuple type with the
    //particles used in a specific context. This one must be
    //specialized for each context
    template<class contextType>
    struct particleTypes {
      static_assert(std::is_base_of<wrapper_context, contextType>::value,
		    "Structure 'particleTypes' used with a 'context'"
		    "not derived from 'wrapper_context'");         
      typedef std::tuple<void> type;
    };
    template<class contextType>
    using particleTypes_t = typename particleTypes<contextType>::type;

    // + Auxiliary structures to handle particle types
    //-------------------------------------------------
    
    // ** Define a structure to get the number of
    //particle types used in a specific context
    template<class contextType>
    struct particleNTypes {
      static_assert(std::is_base_of<wrapper_context, contextType>::value,
		    "Structure 'particleType' used with a 'context'"
		    "not derived from 'wrapper_context'");
      static constexpr const size_t value = std::tuple_size<particleTypes_t<contextType>>::value;
    };
    template<class contextType>
    constexpr const size_t particleNTypes_v = particleNTypes<contextType>::value;
    
    // ** Define struct to get a specific context particle type
    template<class contextType, size_t N>
    struct particleType {
      static_assert(std::is_base_of<wrapper_context, contextType>::value,
		    "Structure 'particleType' used with a 'context'"
		    "not derived from 'wrapper_context'");
      static_assert(N < particleNTypes_v<contextType>,
		    "Structure 'particleType' out of range");
      typedef typename std::tuple_element<N, particleTypes_t<contextType>>::type type;
    };
    template<class contextType, size_t N>
    using particleType_t = typename particleType<contextType, N>::type;

    // ** Define a structure to get a specific context particle state type
    template<class contextType, size_t N>
    struct particleStateType {
      static_assert(std::is_base_of<wrapper_context, contextType>::value,
		    "Structure 'particleStateType' used with a 'context'"
		    "not derived from 'wrapper_context'");
      static_assert(N < particleNTypes_v<contextType>,
		    "Structure 'particleStateType' out of range");      
      typedef typename std::tuple_element<N, particleTypes_t<contextType>>::type::typeState type;
    };
    template<class contextType, size_t N>
    using particleStateType_t = typename particleStateType<contextType, N>::type;
    
    // ** Define a function to search a particle index by type
    template<class contextType, class partType, size_t N = 0>
    constexpr size_t searchType(){
      
      if(N >= particleNTypes<contextType>::value)
	return particleNTypes<contextType>::value;
      
      if(std::is_same<particleType_t<contextType, N>,
	 partType>::value){
	return N;
      }else{
	return searchType<partType, N+1>();
      }
    }
    
    // + Auxiliary structures to handle particle stack types
    //-------------------------------------------------------
    
    // ** Define a function to generate a tuple with stack types for a specific context
    template<class contextType, size_t N = 0>
    constexpr std::enable_if_t<N >= particleNTypes<contextType>::value, std::tuple<void>>
    createStackTuple(){
      return std::tuple<void>();
    }

    template<class contextType, size_t N = 0>
    constexpr auto createStackTuple(std::enable_if_t<N+1 == particleNTypes<contextType>::value, const int> = 0){
      
      typedef particleType_t<contextType, N> ptype;
      //This is the last particle type, return the final stack in a tuple
      return std::tuple<pen_particleStack<typename ptype::typeState>>();      
    }    
    
    template<class contextType, size_t N = 0>
    constexpr auto createStackTuple(std::enable_if_t<N+1 < particleNTypes<contextType>::value, const int> = 0){
      
      typedef particleType_t<contextType, N> ptype;

      //There are more particle types, add the corresponding
      //stack and append the next ones	
      return std::tuple_cat(std::tuple<pen_particleStack<typename ptype::typeState>>(),createStackTuple<contextType, N+1>());      
    }

    // ** Define a structure to get the tuple type with the particle stacks
    //used in a specific context.        
    template<class contextType>
    struct stackTypes {
      static_assert(std::is_base_of<wrapper_context, contextType>::value,
		    "Structure 'stackType' used with a 'context'"
		    "not derived from 'wrapper_context'");
      typedef decltype(createStackTuple<contextType,0>()) type;
    };
    template<class contextType>
    using stackTypes_t = typename stackTypes<contextType>::type;

    // ** Define struct to get a specific context particle type
    template<class contextType, size_t N>
    struct stackType {
      static_assert(std::is_base_of<wrapper_context, contextType>::value,
		    "Structure 'stackType' used with a 'context'"
		    "not derived from 'wrapper_context'");
      static_assert(N < particleNTypes_v<contextType>,
		    "Structure 'stackType' out of range");
      typedef typename std::tuple_element<N, stackTypes_t<contextType>>::type type;
    };
    template<class contextType, size_t N>
    using stackType_t = typename stackType<contextType, N>::type;

    // + Functions to create a callable wrapper with context particles as last arguments
    //-----------------------------------------------------------------------------------

    template <class contextType, class retType, size_t N, class... argTypes>
    inline std::enable_if_t< N == particleNTypes_v<contextType>, std::function<retType(argTypes...)>>
    callableWrapper(const int = 0){
      return std::function<retType(argTypes...)>();
    }
    
    template <class contextType, class retType, size_t N, class... argTypes>
    inline auto callableWrapper(std::enable_if_t< (N != particleNTypes_v<contextType>), const int> = 0){
      
      static_assert(N < particleNTypes_v<contextType>,
		    "callableWrapper: 'N' exceeds the number of available"
		    " particles in this context");
	
      return callableWrapper<contextType, retType, N+1, argTypes..., particleType<contextType, N>>();
    }
    
    // + Structure to instance and construct particles for the specified context
    //---------------------------------------------------------------------------
    
    template<class contextType>
    struct particles{

    public:
      //Store number of particle types
      static constexpr const size_t NTypes =
	particleNTypes_v<contextType>;
      
    private:
      
      // + Functions to register VR in instanced particles
      //--------------------------------------------------------
      
      // Functions to register generic VR
      template < size_t N >
      inline std::enable_if_t<N+1 != NTypes, void>
      _registerGenericVR(const abc_VR<pen_particleState>& vrIn){

	static_assert(N < NTypes,
		      "_registerGenericVR: 'N' exceeds the number of available"
		      " particles in this context. Please, report this error");
	
	getParticle<N>().registerGenericVR(vrIn);
	_registerGenericVR<N+1>(vrIn);
      }

      template < size_t N >
      inline std::enable_if_t<N+1 == NTypes, void> _registerGenericVR(const abc_VR<pen_particleState>& vrIn){
	getParticle<N>().registerGenericVR(vrIn);
      }
      
      // Functions to register specific VR
      template<class stateType, size_t N>
      inline std::enable_if_t<N+1 != NTypes, void>
      _registerSpecificVR(const std::enable_if_t<
			  std::is_same<particleStateType_t<contextType,N>,stateType>::value,
			  abc_VR<stateType>
			  >& vrIn){

	static_assert(N < NTypes,
		      "_registerSpecificVR: 'N' exceeds the number of available"
		      " particles in this context. Please, report this error");
	
	getParticle<N>().registerSpecificVR(vrIn);
	_registerSpecificVR<stateType,N+1>(vrIn);
      }

      template<class stateType, size_t N>
      inline std::enable_if_t<N+1 != NTypes, void>
      _registerSpecificVR(std::enable_if_t<
			  !std::is_same<particleStateType_t<contextType,N>,stateType>::value,
			  const abc_VR<stateType>
			  >& vrIn){

	static_assert(N < NTypes,
		      "_registerSpecificVR: 'N' exceeds the number of available"
		      " particles in this context. Please, report this error");
	
	_registerSpecificVR<stateType,N+1>(vrIn);
      }

      template<class stateType, size_t N>
      inline std::enable_if_t<N+1 == NTypes, void>
      _registerSpecificVR(std::enable_if_t<
			  std::is_same<particleStateType_t<contextType,N>,stateType>::value,
			  const abc_VR<stateType>
			  >& vrIn){
	getParticle<N>().registerSpecificVR(vrIn);
      }

      template<class stateType, size_t N>
      inline std::enable_if_t<N+1 == NTypes, void>
      _registerSpecificVR(std::enable_if_t<
			  !std::is_same<particleStateType_t<contextType,N>,stateType>::value,
			  const abc_VR<stateType>>&){}
      
      // Functions to register an arbitrary number of VR clusters

      //Funciton for non "pen_particleState" state with non empty pack
      template<class stateType, class... otherStateTypes>
      inline std::enable_if_t<sizeof...(otherStateTypes) != 0,void>
      _registerVR(std::enable_if_t<!std::is_same<stateType, pen_particleState>::value, const abc_VR<stateType>>& vr,
		  const abc_VR<otherStateTypes>&... vrRemaining){
	_registerSpecificVR<stateType,0>(vr);
	_registerVR<otherStateTypes...>(vrRemaining...);
      }

      //Function for "pen_particleState" with non empty pack
      template<class stateType, class... otherStateTypes>
      inline std::enable_if_t<sizeof...(otherStateTypes) != 0,void>
      _registerVR(std::enable_if_t<std::is_same<stateType, pen_particleState>::value, const abc_VR<stateType>>& vr,
		  const abc_VR<otherStateTypes>&... vrRemaining){
	_registerGenericVR<0>(vr);
	_registerVR<otherStateTypes...>(vrRemaining...);
      }
      
      //Final funciton for any state with empty pack
      template<class stateType, class... otherStateTypes>
      inline std::enable_if_t<sizeof...(otherStateTypes) == 0,void>
      _registerVR(const std::enable_if_t<!std::is_same<stateType, pen_particleState>::value, const abc_VR<stateType>>& vr,
		  const abc_VR<otherStateTypes>&... vrRemaining){
	_registerSpecificVR<stateType,0>(vr);
      }

      //Function for "pen_particleState" with empty pack
      template<class stateType, class... otherStateTypes>
      inline std::enable_if_t<sizeof...(otherStateTypes) == 0,void>
      _registerVR(const std::enable_if_t<std::is_same<stateType, pen_particleState>::value, const abc_VR<stateType>>& vr,
		  const abc_VR<otherStateTypes>&... vrRemaining){
	_registerGenericVR<0>(vr);
      }
      
      //------------------------------------------------------
      
    protected:


    public:
      
      stackTypes_t<contextType> stackTuple;
      particleTypes_t<contextType> particleTuple;
      
      //Delete default constructor.
      //The constructor must be defined with a specialization
      particles() = delete;
      particles(const contextType&) = delete;

      //Function to get a specific particle
      template <size_t N>
      constexpr particleType_t<contextType, N>& getParticle() {
	return std::get<N>(particleTuple);
      }

      //Function to get a specific particle type
      template <class partType>
      constexpr partType& getParticleType(){
	constexpr const size_t i = searchType<contextType, partType>();

	static_assert(i < NTypes,
		      "getType: Specified particle type not found"
		      "in the context");
	return getParticle<i>();
      }

      //Function to get a specific stack
      template <size_t N>
      constexpr stackType_t<contextType, N>& getStack() {
	return std::get<N>(stackTuple);
      }

      //Function to register generic VR in all particles
      inline void registerGenericVR(const abc_VR<pen_particleState>& vrIn){
	_registerGenericVR<0>();
      }

      //Function to register specific VR in particles with matching state type
      template<class stateType>
      inline void registerSpecificVR(const abc_VR<stateType>& vrIn){
	_registerSpecificVR<stateType, 0>(vrIn);
      }

      //Function to register an arbitrary number of VR clusters
      template<class stateType, class... stateTypes>
      inline void registerVR(const abc_VR<stateType>& vrIn,
			     const abc_VR<stateTypes>&... others){
	_registerVR<stateType, stateTypes...>(vrIn, others...);
      }

      inline void registerVR(){}
      
    };
    
  };

};

#endif
