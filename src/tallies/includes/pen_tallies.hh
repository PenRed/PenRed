
//
//
//    Copyright (C) 2019-2024 Universitat de València - UV
//    Copyright (C) 2019-2024 Universitat Politècnica de València - UPV
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
//        vicent.gimenez.alventosa@gmail.com (Vicent Giménez Alventosa)
//        vicente.gimenez@uv.es (Vicente Giménez Gómez)
//    
//


#ifndef __PEN_TALLIES__
#define __PEN_TALLIES__

#include <cstdio>
#include <string>
#include <vector>
#include <algorithm>
#include <thread>
#include <stdexcept>
#include <functional>
#include <tuple>
#include <type_traits>
#include <typeinfo>
#include <utility>

#include "pen_auxiliar.hh"

// ******************************* MPI ************************************ //
#ifdef _PEN_USE_MPI_
#include "mpi.h"
#endif
// ***************************** MPI END ********************************** //


#include "pen_classes.hh"
#include "pen_math.hh"
#include "instantiator.hh"
#include "pen_data.hh"
#include "pen_dump.hh"
#include "pen_image.hh"
#include "pen_geometries.hh"
#include "pen_states.hh"


template <class stateType> class pen_genericTally;

// Auxiliar structures
//////////////////////////

struct tally_StepData{
  double dsef;
  double dstot;
  double softDE;
  double softX;
  double softY;
  double softZ;
  unsigned originIBODY;
  unsigned originMAT;

  template <class particleType>
  inline void update(const particleType& particle){
    dsef = particle.DSef();
    dstot = particle.DStot();
    originMAT = particle.lastMat();
    originIBODY = particle.lastBody();
  }  
  
};

namespace penred{
  namespace tally{

    template <size_t N>
    struct Dim{
      static constexpr size_t value = N;
    };

    // ++ Results types extraction
    // Primary template (handles empty case)
    template <typename... Args>
    struct ResultsTypeSelector {
      using type = void;  // Default to void
      static_assert(sizeof...(Args) == 0,
		    "Tally return types require either no types or std::pair<T,dim>");
    };
    
    // Specialization for valid pairs
    template <typename T, typename D>
    struct ResultsTypeSelector<std::pair<T, D>> {
      static constexpr size_t dim = D::value;
      using value_type = T;
      
      using type = std::conditional_t<
        (std::is_arithmetic<T>::value && dim > 0),
        measurements::results<T, dim>,
        std::vector<T>
	>;
    };

    // + Results Tuple builder
    template <typename... Pairs>
    struct ResultsTupleBuilder {
      using type = std::tuple<typename ResultsTypeSelector<Pairs>::type...>;
    };

    // Empty case specialization
    template <>
    struct ResultsTupleBuilder<> {
      using type = std::tuple<>;
    };    
        
    template <typename... Args>
    using ResultsType = typename ResultsTupleBuilder<Args...>::type;
    
    // + Results generation tuple builder
    template <typename... Pairs>
    struct ResultsGeneratorTupleBuilder{
      using type = std::tuple<std::function<typename ResultsTypeSelector<Pairs>::type (const unsigned long long)>...>;
    };

    // Empty case specialization
    template <>
    struct ResultsGeneratorTupleBuilder<> {
      using type = std::tuple<>;
    };    

    template <typename... Args>
    using ResultsGeneratorType = typename ResultsGeneratorTupleBuilder<Args...>::type;
    
    // ++ Printing results types functions

    template <typename T>
    constexpr const char* typeNames(){
      if(std::is_same<T, double>::value)
	return "double";
      if(std::is_same<T, float>::value)
	return "float";
      if(std::is_same<T, int>::value)
	return "int";
      if(std::is_same<T, unsigned>::value)
	return "unsigned";
      if(std::is_same<T, long int>::value)
	return "long int";
      if(std::is_same<T, long unsigned>::value)
	return "long unsigned";
      if(std::is_same<T, long long int>::value)
	return "long long int";
      if(std::is_same<T, long long unsigned>::value)
	return "long long unsigned";
      if(std::is_same<T, char>::value)
	return "char";
      if(std::is_same<T, unsigned char>::value)
	return "unsigned char";
      return nullptr;
    }

    template <typename T>
    typename std::enable_if_t<measurements::is_results<T>::value, void>
    resultsTypeInfo(){
      using ValueType = typename measurements::is_results<T>::value_type;
      constexpr const char* mapName = typeNames<ValueType>();
      if(mapName != nullptr){
	std::cout << "   - measurements::results<" 
		  << mapName << ", " 
		  << measurements::is_results<T>::dimension << ">\n";
      }
      else{
	std::cout << "   - measurements::results<" 
		  << typeid(ValueType).name() << ", " 
		  << measurements::is_results<T>::dimension << ">\n";
      }
    }

    template <typename T>
    typename std::enable_if_t<measurements::is_vector<T>::value, void>
    resultsTypeInfo(){
      constexpr const char* mapName = typeNames<measurements::vector_value_type_t<T>>();
      if(mapName != nullptr){
	std::cout << "   - std::vector<" 
		  << mapName << ">\n";
      }
      else{
	std::cout << "   - std::vector<" 
		  << typeid(T).name() << ">\n";
      }
    }    
    
    template <typename T>
    typename std::enable_if_t<!measurements::is_results<T>::value && !measurements::is_vector<T>::value, void>
    resultsTypeInfo(){
      std::cout << "   - " << typeid(T).name() << "\n";	
    }

    // Recursive results type printer
    template <size_t I = 0, typename Tuple>
    typename std::enable_if<I == std::tuple_size<Tuple>::value>::type
    printResultsTypes(const Tuple&) {
      if(I == 0){
	std::cout << "   - Null \n" << std::endl;
      }
    }
    
    template <size_t I = 0, typename Tuple>
    typename std::enable_if<I < std::tuple_size<Tuple>::value>::type
    printResultsTypes(const Tuple& t) {
      //Print type information
      resultsTypeInfo<std::tuple_element_t<I, Tuple>>();
      // Process next element
      printResultsTypes<I+1>(t);
    }

    // Print results types for a specific tally
    template<class TallyType>
    void printResultsTypes(){
      std::cout << TallyType::___ID << ":\n";
      const typename TallyType::ResultsTypes aux;
      printResultsTypes(aux);
    }
    
  } //End of tally namespace
} //End of penred namespace

// Register macros
///////////////////////

#define DECLARE_TALLY(Class, State, ID, ...)				\
  public:								\
  using ResultsTypes = penred::tally::ResultsType<__VA_ARGS__>;		\
  using ResultsMapType = std::map<std::string, ResultsTypes>;		\
  using ResultsGeneratorType = penred::tally::ResultsGeneratorType<__VA_ARGS__>; \
  ResultsGeneratorType resultsGenerators;				\
  static int registerStatus();						\
  const char* readID() const;						\
  static constexpr const char* ___ID = static_cast<const char *>(#ID);	\
  static constexpr const char* tallyID() { return ___ID; }		\
  static volatile const int ___register_return;				\
  inline int sum(const pen_genericTally<State>& sumtally){		\
    const Class& derived = dynamic_cast<const Class&>(sumtally);	\
    return sumTally(derived);						\
  }									\
  inline int shareConfig(const pen_genericTally<State>& sharingTally){	\
    const Class& derived = dynamic_cast<const Class&>(sharingTally);	\
    return sharedConfig(derived);					\
  }									\
  static inline void printResultsTypes(){				\
    penred::tally::printResultsTypes<Class>();				\
  }									\
  /* Results generators */						\
  template<size_t I>							\
  std::tuple_element_t<I, ResultsTypes> generateResults(const unsigned long long __nhists){ \
  auto generator = std::get<I>(resultsGenerators);			\
  if(generator)								\
    return generator(__nhists);							\
  else									\
    return std::tuple_element_t<I, ResultsTypes>();			\
  }									\
   template<size_t I>							\
   inline void updateResult(ResultsTypes& r, const unsigned long long __nhists){ \
     std::get<I>(r) = generateResults<I>(__nhists);				\
   }									\
   template <size_t... I>						\
   inline void updateResults(ResultsTypes& r, const unsigned long long __nhists, const std::index_sequence<I...>&){ \
     int dummy[] = {0, (updateResult<I>(r, __nhists), 0)...};		\
     (void)dummy;							\
     (void)__nhists;							\
   }									\
   inline void updateResults(ResultsTypes& r, const unsigned long long __nhists){	\
     updateResults(r, __nhists, std::make_index_sequence<std::tuple_size<ResultsTypes>::value>{}); \
   }									\
   inline ResultsTypes getResults(const unsigned long long __nhists){	\
     ResultsTypes r;							\
     updateResults(r, __nhists);							\
     return r;								\
   }									\
   template<size_t I>							\
   inline void setResultsGenerator(std::tuple_element_t<I, ResultsGeneratorType> f){ \
     std::get<I>(resultsGenerators) = f;				\
   }									\
private:

#define REGISTER_COMMON_TALLY(Class) \
  volatile const int Class::___register_return = pen_commonTallyCluster::addTally<Class>(Class::tallyID()); \
  int Class::registerStatus() { return ___register_return;}		\
  const char* Class::readID() const { return tallyID();}

// Enumeration flags
//////////////////////////

enum __usedFunc{
	      USE_BEGINSIM     = 1 << 0,
	      USE_ENDSIM       = 1 << 1,
	      USE_SAMPLEDPART  = 1 << 2,
	      USE_ENDHIST      = 1 << 3,
	      USE_MOVE2GEO     = 1 << 4,
	      USE_BEGINPART    = 1 << 5,
	      USE_ENDPART      = 1 << 6,
	      USE_JUMP         = 1 << 7,
	      USE_STEP         = 1 << 8,
	      USE_INTERFCROSS  = 1 << 9,
	      USE_MATCHANGE    = 1 << 10,
	      USE_KNOCK        = 1 << 11,
	      USE_LOCALEDEP    = 1 << 12,
	      USE_LASTHIST     = 1 << 13,
	      INVALID_FUNCTION = 1 << 14
		
};


inline __usedFunc operator|(__usedFunc a, __usedFunc b)
{return static_cast<__usedFunc>(static_cast<int>(a) | static_cast<int>(b));}

// Tally class
//////////////////////////

template <class stateType>
class pen_genericTally : public penred::logs::logger{

private:

  unsigned nthread;
  std::string name;
  std::string OutputDirPath;
  const __usedFunc usedFunctions;

  //Create a vector of image exporters
  std::vector<pen_imageExporter> imageExporters;

  //Create an array to store particle stacks
  std::array<const abc_particleStack*, constants::nParTypes> stacks;
protected:

  //Create dump instance
  pen_dump dump;

  //Include sub dump
  inline void addSubDump(pen_genericTally<stateType>& subTally){
    dump.toDump(subTally.dump);
  }

  //Include measurement in dump
  template<class valType, size_t dim>
  inline void toDump(penred::measurements::measurement<valType, dim>& m){
    dump.toDump(m.getData().data(), m.getNBins());
    dump.toDump(m.getData2().data(), m.getNBins());
  }

  //Create filenames with output dir path, thread and MPI number
  std::string createFilename(const char* filename) const{

    // ******************************* MPI ********************************** //
#ifdef _PEN_USE_MPI_
    int rank;
    //Add the MPI rank number to filename
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    return std::string(OutputDirPath + name + std::string("-MPI") +
		       std::to_string(rank) + std::string("-th") +
		       std::to_string(nthread) + std::string("-") +
		       std::string(filename));

    // ***************************** MPI END ******************************** //
#else
    
    return std::string(OutputDirPath + name + std::string("-th") +
		       std::to_string(nthread) + std::string("-") +
		       std::string(filename));

#endif    
  }  
  
  //Create private fopen method to append tally name and thread
  //identifier to generated and read files
  inline FILE* fopen(const char* filename, const char* mode) const{

    if(filename == nullptr || mode == nullptr)
      return nullptr;

    std::string finalFilename;
    FILE* fout = nullptr;
    
  // ******************************* MPI ************************************ //
#ifdef _PEN_USE_MPI_
    int rank;
    //Add the MPI rank number to filename
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    finalFilename = OutputDirPath + name + std::string("-MPI") +
      std::to_string(rank) + std::string("-th") + std::to_string(nthread) +
      std::string("-") + std::string(filename);

  // ***************************** MPI END ********************************** //
#else
    
    finalFilename = OutputDirPath + name + std::string("-th") +
      std::to_string(nthread) + std::string("-") + std::string(filename);

#endif


    fout = ::fopen(finalFilename.c_str(),mode);
    //Check if the file can be opened
    if(fout == nullptr){
      //Try to open the same file name but with a 'not-found' prefix
      printf("pen_genericTally: fopen: Error: Unable to open "
	     "file ('%s'): \n%s\n",
	     mode, finalFilename.c_str());
      
  // ******************************* MPI ************************************ //
#ifdef _PEN_USE_MPI_

      finalFilename = std::string("not-found-") + name + std::string("-MPI") +
	std::to_string(rank) + std::string("-th") + std::to_string(nthread) +
	std::string("-") + std::string(filename);
      
  // ***************************** MPI END ********************************** //
#else
      
      finalFilename = std::string("not-found-") + name + std::string("-th") +
	std::to_string(nthread) + std::string("-") + std::string(filename);
      
#endif

      printf("Trying to redirect to: '%s'\n",finalFilename.c_str());
      fflush(stdout);

      fout = ::fopen(finalFilename.c_str(),mode);
    }
    //Return resulting file pointer
    return fout;
  }

public:

  pen_genericTally(__usedFunc functionsToUse) : nthread(0),
						name(""),
						usedFunctions(functionsToUse)
  {
    std::fill(stacks.begin(), stacks.end(), nullptr);
  }

  virtual const char* readID() const = 0;
  inline const std::string& readName() const { return name;}
  inline unsigned getThread() const {return nthread;}
  inline const std::string& readOutputDirPath() const { return OutputDirPath;}
  inline void setName(const char* newName) { name.assign(newName);}
  inline void setName(const std::string& newName) { name.assign(newName);}
  inline void setOutputDirPath(const char* newName) { OutputDirPath.assign(newName);}
  inline void setOutputDirPath(const std::string& newName) { OutputDirPath.assign(newName);}
  inline void setThread(const unsigned threadNum){nthread = threadNum;}
  inline void setStack(const pen_KPAR kpar, const abc_particleStack* stack){stacks[kpar] = stack;}
  inline const abc_particleStack* readStack(const pen_KPAR kpar) const {return stacks[kpar];}

  template<class T>
  inline void addImage(const char* filename,
		       const unsigned nDim,
		       const unsigned* elements,
		       const float* delements,
		       const double* origin,
		       std::function<T(unsigned long long, size_t)> fin){
    imageExporters.emplace(imageExporters.end(), fin);
    imageExporters.back().setDimensions(nDim,elements,delements);
    imageExporters.back().setOrigin(origin);
    imageExporters.back().baseName = createFilename(filename);
  }

  template<class T>
  inline void addImage(const char* filename,
		       const unsigned nDim,
		       const unsigned* elements,
		       const float* delements,
		       const double* origin,
		       std::function<T(unsigned long long, size_t, T&)> fin){
    imageExporters.emplace(imageExporters.end(), fin);
    imageExporters.back().setDimensions(nDim,elements,delements);
    imageExporters.back().setOrigin(origin);
    imageExporters.back().baseName = createFilename(filename);
  }  

  inline void exportImage(const unsigned long long nhists,
			  const pen_imageExporter::formatTypes format) const {
    for(auto& exporter : imageExporters)
      exporter.exportImage(nhists,format);
  }
  
  void getUsedMethods(std::vector<std::string>& funcNames) const{
    if(has_beginSim())
      funcNames.push_back(std::string("beginSim"));
    if(has_endSim())
      funcNames.push_back(std::string("endSim"));  
    if(has_sampledPart())
      funcNames.push_back(std::string("sampledPart"));  
    if(has_endHist())
      funcNames.push_back(std::string("endHist"));
    if(has_move2geo())
      funcNames.push_back(std::string("move2geo"));
    if(has_beginPart())
      funcNames.push_back(std::string("beginPart"));  
    if(has_endPart())
      funcNames.push_back(std::string("endPart"));  
    if(has_jump())
      funcNames.push_back(std::string("jump"));      
    if(has_step())
      funcNames.push_back(std::string("step"));
    if(has_interfCross())
      funcNames.push_back(std::string("interfCross"));      
    if(has_matChange())
      funcNames.push_back(std::string("matChange"));      
    if(has_knock())
      funcNames.push_back(std::string("knock"));      
    if(has_localEdep())
      funcNames.push_back(std::string("localEdep"));
    if(has_lastHist())
      funcNames.push_back(std::string("lastHist"));    
  }
  
  virtual void tally_beginSim(){
    char error[1000];
    sprintf(error,"Error on tally %s of type %s: \n Trying to call 'tally_beginSim' but has not been implemented.\n",readName().c_str(),readID());
    printf("%s\n",error);
    fflush(stdout);
    throw std::bad_function_call();
  }
  virtual void tally_endSim(const unsigned long long /*nhist*/){
    char error[1000];
    sprintf(error,"Error on tally %s of type %s: \n Trying to call 'tally_endSim' but has not been implemented.\n",readName().c_str(),readID());
    printf("%s\n",error);
    fflush(stdout);
    throw std::bad_function_call();
  }
  
  virtual void tally_sampledPart(const unsigned long long /*nhist*/,
				 const unsigned long long /*dhist*/,
				 const unsigned /*kdet*/,
				 const pen_KPAR /*kpar*/,
				 const stateType& /*state*/){
    char error[1000];
    sprintf(error,"Error on tally %s of type %s: \n Trying to call 'tally_sampledPart' but has not been implemented.\n",readName().c_str(),readID());
    printf("%s\n",error);
    fflush(stdout);
    throw std::bad_function_call();
  }
  
  virtual void tally_endHist(const unsigned long long /*nhist*/){
    char error[1000];
    sprintf(error,"Error on tally %s of type %s: \n Trying to call 'tally_endHist' but has not been implemented.\n",readName().c_str(),readID());
    printf("%s\n",error);
    fflush(stdout);
    throw std::bad_function_call();
  }

  virtual void tally_move2geo(const unsigned long long /*nhist*/,
			      const unsigned /*kdet*/,
			      const pen_KPAR /*kpar*/,
			      const stateType& /*state*/,
			      const double /*dsef*/,
			      const double /*dstot*/){
    char error[1000];
    sprintf(error,"Error on tally %s of type %s: \n Trying to call 'tally_move2geo' but has not been implemented.\n",readName().c_str(),readID());
    printf("%s\n",error);
    fflush(stdout);
    throw std::bad_function_call();
  }
  
  virtual void tally_beginPart(const unsigned long long /*nhist*/,
			       const unsigned /*kdet*/,
			       const pen_KPAR /*kpar*/,
			       const stateType& /*state*/){
    char error[1000];
    sprintf(error,"Error on tally %s of type %s: \n Trying to call 'tally_beginPart' but has not been implemented.\n",readName().c_str(),readID());
    printf("%s\n",error);
    fflush(stdout);
    throw std::bad_function_call();
  }
  
  virtual void tally_endPart(const unsigned long long /*nhist*/,
			     const pen_KPAR /*kpar*/,
			     const stateType& /*state*/){
    char error[1000];
    sprintf(error,"Error on tally %s of type %s: \n Trying to call 'tally_endPart' but has not been implemented.\n",readName().c_str(),readID());
    printf("%s\n",error);
    fflush(stdout);
    throw std::bad_function_call();
  }
  
  virtual void tally_localEdep(const unsigned long long /*nhist*/,
			       const pen_KPAR /*kpar*/,
			       const stateType& /*state*/,
			       const double /*dE*/){
    char error[1000];
    sprintf(error,"Error on tally %s of type %s: \n Trying to call 'tally_localEdep' but has not been implemented.\n",readName().c_str(),readID());
    printf("%s\n",error);
    fflush(stdout);
    throw std::bad_function_call();
  }
  
  virtual void tally_step(const unsigned long long /*nhist*/,
			  const pen_KPAR /*kpar*/,
			  const stateType& /*state*/,
			  const tally_StepData& /*stepData*/){
    char error[1000];
    sprintf(error,"Error on tally %s of type %s: \n Trying to call 'tally_step' but has not been implemented.\n",readName().c_str(),readID());
    printf("%s\n",error);
    fflush(stdout);
    throw std::bad_function_call();
  }
  
  virtual void tally_interfCross(const unsigned long long /*nhist*/,
				 const unsigned /*kdet*/,
				 const pen_KPAR /*kpar*/,
				 const stateType& /*state*/){
    char error[1000];
    sprintf(error,"Error on tally %s of type %s: \n Trying to call 'tally_interfCross' but has not been implemented.\n",readName().c_str(),readID());
    printf("%s\n",error);
    fflush(stdout);
    throw std::bad_function_call();
  }
  
  virtual void tally_matChange(const unsigned long long /*nhist*/,
			       const pen_KPAR /*kpar*/,
			       const stateType& /*state*/,
			       const unsigned /*prevMat*/){
    char error[1000];
    sprintf(error,"Error on tally %s of type %s: \n Trying to call 'tally_matChange' but has not been implemented.\n",readName().c_str(),readID());
    printf("%s\n",error);
    fflush(stdout);
    throw std::bad_function_call();
  }
  
  virtual void tally_jump(const unsigned long long /*nhist*/,
			  const pen_KPAR /*kpar*/,
			  const stateType& /*state*/,
			  const double /*ds*/){
    char error[1000];
    sprintf(error,"Error on tally %s of type %s: \n Trying to call 'tally_jump' but has not been implemented.\n",readName().c_str(),readID());
    printf("%s\n",error);
    fflush(stdout);
    throw std::bad_function_call();
  }
  
  virtual void tally_knock(const unsigned long long /*nhist*/,
			   const pen_KPAR /*kpar*/,
			   const stateType& /*state*/,
			   const int /*icol*/){
    char error[1000];
    sprintf(error,"Error on tally %s of type %s: \n Trying to call 'tally_knock' but has not been implemented.\n",readName().c_str(),readID());
    printf("%s\n",error);
    fflush(stdout);
    throw std::bad_function_call();
  }

  virtual void tally_lastHist(const unsigned long long /*lasthist*/){
    char error[1000];
    sprintf(error,"Error on tally %s of type %s: \n Trying to call 'tally_lastHist' but has not been implemented.\n",readName().c_str(),readID());
    printf("%s\n",error);
    fflush(stdout);
    throw std::bad_function_call();
  }
  
  virtual int configure(const wrapper_geometry& /*geometry*/,
			const abc_material* const /*materials*/[constants::MAXMAT],
			const pen_parserSection& /*config*/,
			const unsigned /*verbose*/) = 0;
  
  virtual void saveData(const unsigned long long /*nhist*/) const = 0;
  virtual int sum(const pen_genericTally<stateType>&) = 0;
  virtual int shareConfig(const pen_genericTally<stateType>&) = 0;
  int sharedConfig(const pen_genericTally<stateType>&){return 0;}
  virtual void flush() = 0;
  
  inline __usedFunc readFlags() const {return usedFunctions;}
  
  inline bool has_beginSim() const{
    return ((usedFunctions & USE_BEGINSIM) == USE_BEGINSIM);
  }
  inline bool has_endSim() const{
    return ((usedFunctions & USE_ENDSIM) == USE_ENDSIM);
  }
  inline bool has_sampledPart() const{
    return ((usedFunctions & USE_SAMPLEDPART) == USE_SAMPLEDPART);
  }
  inline bool has_endHist() const{
    return ((usedFunctions & USE_ENDHIST) == USE_ENDHIST);
  }
  inline bool has_move2geo() const{
    return ((usedFunctions & USE_MOVE2GEO) == USE_MOVE2GEO);
  }
  inline bool has_beginPart() const{
    return ((usedFunctions & USE_BEGINPART) == USE_BEGINPART);
  }
  inline bool has_endPart() const{
    return ((usedFunctions & USE_ENDPART) == USE_ENDPART);
  }
  inline bool has_jump() const{
    return ((usedFunctions & USE_JUMP) == USE_JUMP);
  }
  inline bool has_step() const{
    return ((usedFunctions & USE_STEP) == USE_STEP);
  }
  inline bool has_interfCross() const{
    return ((usedFunctions & USE_INTERFCROSS) == USE_INTERFCROSS);
  }
  inline bool has_matChange() const{
    return ((usedFunctions & USE_MATCHANGE) == USE_MATCHANGE);
  }
  inline bool has_knock() const{
    return ((usedFunctions & USE_KNOCK) == USE_KNOCK);
  }
  inline bool has_localEdep() const{
    return ((usedFunctions & USE_LOCALEDEP) == USE_LOCALEDEP);
  }
  inline bool has_lastHist() const{
    return ((usedFunctions & USE_LASTHIST) == USE_LASTHIST);
  }
  

  inline int writeDump(unsigned char*& pout, size_t& pos, const unsigned verbose = 0) {
    //Flush tally
    flush();
    //Dump it
    return dump.dump(pout,pos,0,verbose);
  }

  inline int readDump(const unsigned char* const pin, size_t& pos, const unsigned verbose = 0){
    return dump.read(pin,pos,verbose);
  }
  
  virtual ~pen_genericTally(){}
  
};

// Auxiliar functions
///////////////////////

bool __tallySort (pen_genericTally<pen_particleState>* i,
		  pen_genericTally<pen_particleState>* j);

namespace penred{
  namespace tally{
    
    //Check registered types
    template <class stateType>
    bool checkRegistered(const unsigned verbose);

    // ++ Downcasting via unique ID check
    template<class TallyType>
    typename std::enable_if_t<std::is_base_of<pen_genericTally<pen_particleState>, TallyType>::value, TallyType>*
    downcast(pen_genericTally<pen_particleState>* ptally){
      if(ptally->readID() == TallyType::tallyID())
	return static_cast<TallyType*>(ptally);
      else
	return nullptr;
    }        
  }
}

// Tally cluster classes
//////////////////////////

class pen_commonTallyCluster : public penred::logs::logger{

private:
  unsigned nthread;
protected:

  typedef pen_genericTally<pen_particleState>* ptallyType;
  typedef std::vector<ptallyType>::iterator tallyIterator;
  
  //Tallies instantiators
  static instantiator<pen_genericTally<pen_particleState>>& genericTallies();
  
  //Registered tallies vector
  std::vector<pen_genericTally<pen_particleState>*> tallies;

  //Vectors containing tally pointers with specific collection functions
  std::vector<pen_genericTally<pen_particleState>*> tallies_beginSim;
  std::vector<pen_genericTally<pen_particleState>*> tallies_endSim;
  std::vector<pen_genericTally<pen_particleState>*> tallies_sampledPart;
  std::vector<pen_genericTally<pen_particleState>*> tallies_endHist;
  std::vector<pen_genericTally<pen_particleState>*> tallies_move2geo;
  std::vector<pen_genericTally<pen_particleState>*> tallies_beginPart;
  std::vector<pen_genericTally<pen_particleState>*> tallies_endPart;
  std::vector<pen_genericTally<pen_particleState>*> tallies_localEdep;
  std::vector<pen_genericTally<pen_particleState>*> tallies_step;
  std::vector<pen_genericTally<pen_particleState>*> tallies_interfCross;
  std::vector<pen_genericTally<pen_particleState>*> tallies_matChange;
  std::vector<pen_genericTally<pen_particleState>*> tallies_jump;
  std::vector<pen_genericTally<pen_particleState>*> tallies_knock;
  std::vector<pen_genericTally<pen_particleState>*> tallies_lastHist;

  int configStatus;

  pen_commonTallyCluster* mpiBuffer;

  int createTally(const char* OutDir, const char* ID,
		  const char* name,
		  const wrapper_geometry& geometry,
		  const abc_material* const materials[constants::MAXMAT], 
		  const pen_parserSection& config,
		  const unsigned verbose);  
public:

  std::string name;

  pen_commonTallyCluster() : nthread(0),
			     configStatus(0),
			     mpiBuffer(nullptr),
			     name("unnamed")
  {}

  inline int configureStatus(){return configStatus;}
  
  inline unsigned numTallies(){return tallies.size();}
  
  static std::string typesList(){
    return genericTallies().typesList();
  }

  static std::string typesList(std::vector<std::string>& list){
    return genericTallies().typesList(list);
  }

  static pen_genericTally<pen_particleState>* createInstance(const char* type){
    return genericTallies().createInstance(type);
  }  
  static pen_genericTally<pen_particleState>* createInstance(const std::string type){
    return genericTallies().createInstance(type.c_str());
  }

  inline unsigned getThread() const {return nthread;}
  
  template <class subclass>
  static int addTally(const char* typeID){
    return genericTallies().addSubType<subclass>(typeID);
  }

  inline int getTallyIndex(const std::string& tallyName){
    for(unsigned i = 0; i < tallies.size(); ++i){
      if(tallies[i]->readName().compare(tallyName) == 0)
	return static_cast<int>(i);
    }
    return -1;
  }
  inline std::string getTallyName(const size_t i){
    if(i >= tallies.size())
      return std::string();

    return tallies[i]->readName();
  }

  //Define a set of functions to get tally results
  template<class TallyType>
  inline bool isCreated(const std::string& tallyName){
    for(pen_genericTally<pen_particleState>* t : tallies){
      TallyType* derived = penred::tally::downcast<TallyType>(t);
      if(derived != nullptr){
	if(derived->readName().compare(tallyName) == 0){
	  return true;
	}
      }
    }
    //Not found
    return false;
  }
  template<class TallyType>
  inline bool isCreated(const size_t i){
    if(i >= tallies.size())
      return false;
    
    TallyType* derived = penred::tally::downcast<TallyType>(tallies[i]);
    if(derived != nullptr){
	return true;
    }
    
    //Not found
    return false;
  }
  
  template<class TallyType>
  inline typename TallyType::ResultsTypes getResults(const std::string& tallyName, const unsigned long long nhists){
    for(pen_genericTally<pen_particleState>* t : tallies){
      TallyType* derived = penred::tally::downcast<TallyType>(t);
      if(derived != nullptr){
	if(derived->readName().compare(tallyName) == 0){
	  return derived->getResults(nhists);
	}
      }
    }
    //Not found
    return typename TallyType::ResultsTypes();
  }
  template<class TallyType>
  inline typename TallyType::ResultsTypes getResults(const size_t i, const unsigned long long nhists){

    if(i >= tallies.size())
      return typename TallyType::ResultsTypes();

    TallyType* derived = penred::tally::downcast<TallyType>(tallies[i]);
    if(derived != nullptr){
      return derived->getResults(nhists);
    }

    //Not found
    return typename TallyType::ResultsTypes();
  }

  std::thread configure_async(const wrapper_geometry* geometry,
			      const abc_material* const materials[constants::MAXMAT],
			      const unsigned threadNum,
			      const pen_parserSection& config,
			      const unsigned verbose = 0);
  
  void configure(const wrapper_geometry* geometry,
		 const abc_material* const materials[constants::MAXMAT],
		 const unsigned threadNum,
		 const pen_parserSection config,
		 const unsigned verbose = 0);

  inline void setStack(const pen_KPAR kpar, const abc_particleStack* stack){
    for(tallyIterator i = tallies.begin();
	i != tallies.end(); ++i)
      (*i)->setStack(kpar,stack);
  }
  
  inline void run_beginSim(){
    for(tallyIterator i = tallies_beginSim.begin();
	i != tallies_beginSim.end(); ++i)
      (*i)->tally_beginSim();
  }

  inline void run_endSim(const unsigned long long nhist){
    for(tallyIterator i = tallies_endSim.begin();
	i != tallies_endSim.end(); ++i)
      (*i)->tally_endSim(nhist);
  }
  
  inline void run_sampledPart(const unsigned long long nhist,
			      const unsigned long long dhist,
			      const unsigned kdet,
			      const pen_KPAR kpar,
			      const pen_particleState& state){

    for(tallyIterator i = tallies_sampledPart.begin();
	i != tallies_sampledPart.end(); ++i)
      (*i)->tally_sampledPart(nhist,dhist,kdet,kpar,state);
  }
  
  inline void run_endHist(const unsigned long long nhist){

    for(tallyIterator i = tallies_endHist.begin();
	i != tallies_endHist.end(); ++i)
      (*i)->tally_endHist(nhist);
  }

  inline void run_move2geo(const unsigned long long nhist,
			   const unsigned kdet,
			   const pen_KPAR kpar,
			   const pen_particleState& state,
			   const double dsef,
			   const double dstot){
    
    for(tallyIterator i = tallies_move2geo.begin();
	i != tallies_move2geo.end(); ++i)
      (*i)->tally_move2geo(nhist,kdet,kpar,state,dsef,dstot);
  }
  
  inline void run_beginPart(const unsigned long long nhist,
			    const unsigned kdet,
			    const pen_KPAR kpar,
			    const pen_particleState& state){

    for(tallyIterator i = tallies_beginPart.begin();
	i != tallies_beginPart.end(); ++i)
      (*i)->tally_beginPart(nhist,kdet,kpar,state);
  }
  
  inline void run_endPart(const unsigned long long nhist,
			  const pen_KPAR kpar,
			  const pen_particleState& state){

    for(tallyIterator i = tallies_endPart.begin();
	i != tallies_endPart.end(); ++i)
      (*i)->tally_endPart(nhist,kpar,state);
  }
  
  inline void run_localEdep(const unsigned long long nhist,
			const pen_KPAR kpar,
			const pen_particleState& state,
			const double dE){

    for(tallyIterator i = tallies_localEdep.begin();
	i != tallies_localEdep.end(); ++i)
      (*i)->tally_localEdep(nhist,kpar,state,dE);
  }
  
  inline void run_step(const unsigned long long nhist,
		       const pen_KPAR kpar,
		       const pen_particleState& state,
		       const tally_StepData& stepData){

    for(tallyIterator i = tallies_step.begin();
	i != tallies_step.end(); ++i)
      (*i)->tally_step(nhist,kpar,state,stepData);    
  }
  
  inline void run_interfCross(const unsigned long long nhist,
			      const unsigned kdet,
			      const pen_KPAR kpar,
			      const pen_particleState& state){

    for(tallyIterator i = tallies_interfCross.begin();
	i != tallies_interfCross.end(); ++i)
      (*i)->tally_interfCross(nhist,kdet,kpar,state);
  }
  
  inline void run_matChange(const unsigned long long nhist,
			    const pen_KPAR kpar,
			    const pen_particleState& state,
			    const unsigned prevMat){

    for(tallyIterator i = tallies_matChange.begin();
	i != tallies_matChange.end(); ++i)
      (*i)->tally_matChange(nhist,kpar,state,prevMat);
  }
  
  inline void run_jump(const unsigned long long nhist,
		       const pen_KPAR kpar,
		       const pen_particleState& state,
		       const double ds){

    for(tallyIterator i = tallies_jump.begin();
	i != tallies_jump.end(); ++i)
      (*i)->tally_jump(nhist,kpar,state,ds);
  }
  
  inline void run_knock(const unsigned long long nhist,
			const pen_KPAR kpar,
			const pen_particleState& state,
			const int icol){

    for(tallyIterator i = tallies_knock.begin();
	i != tallies_knock.end(); ++i)
      (*i)->tally_knock(nhist,kpar,state,icol);
  }

  inline void run_lastHist(const unsigned long long lastHist){

    for(tallyIterator i = tallies_lastHist.begin();
	i != tallies_lastHist.end(); ++i)
      (*i)->tally_lastHist(lastHist);
  }
  
  inline void saveData(const unsigned long long nhist, const bool doflush = true){

    for(tallyIterator i = tallies.begin();
	i != tallies.end(); ++i){
      if(doflush)
	(*i)->flush();
      (*i)->saveData(nhist);
    }
  }

  inline void exportImage(const unsigned long long nhist,
			  const pen_imageExporter::formatTypes format,
			  const bool doflush = true){

    for(tallyIterator i = tallies.begin();
	i != tallies.end(); ++i){
      if(doflush)
	(*i)->flush();
      (*i)->exportImage(nhist,format);
    }
  }  
  
  int writeDump(unsigned char*& pdump,
		size_t& dim,
		const unsigned long long nhist,
		const int seed1, const int seed2,
		const int lastSource,
		const unsigned long long sourceHists,
		const unsigned verbose = 0);
  int readDump(const unsigned char* const pdump,
	       size_t& pos,
	       unsigned long long& nhist,
	       int& seed1, int& seed2,
	       int& lastSource,
	       unsigned long long& sourceHists,
	       const unsigned verbose = 0);

  int dump2file(const char* filename,
		const unsigned long long nhist,
		const int seed1, const int seed2,
		const int lastSource,
		const unsigned long long sourceHists,
		const unsigned verbose = 0);

  int readDumpfile(const char* filename,
		   unsigned long long& nhist,
		   int& seed1, int& seed2,
		   int& lastSource,
		   unsigned long long& sourceHists,
		   const unsigned verbose = 0);
  
  int sum(pen_commonTallyCluster& cluster, const unsigned verbose = 0);
  int shareConfig(const pen_commonTallyCluster& cluster, const unsigned verbose = 0);
  
  // ******************************* MPI ************************************ //
#ifdef _PEN_USE_MPI_
  int reduceMPI(unsigned long long& simulatedHists,
		const MPI_Comm comm,
		const unsigned verbose);
#endif
  // ***************************** MPI END ********************************** //
  
  void clear();

  ~pen_commonTallyCluster();
  
};

//Include defined tallies
#include "genericTallies.hh"
#include "specificTallies.hh"

// Define a complete results map handler class

namespace penred{
  namespace tally{

    // ++ Results type printing
    template<size_t T>
    typename std::enable_if_t<T >= std::tuple_size<typesGenericTallies>::value, void>
    printTalliesResultsTypes(){}
    
    template<size_t T = 0>
    typename std::enable_if_t<T < std::tuple_size<typesGenericTallies>::value, void>
    printTalliesResultsTypes(){
      using tupleType = typename std::tuple_element<T, typesGenericTallies>::type;
      
      //Print information for this tally type
      printResultsTypes<tupleType>();
      
      //Go to the next tally type
      printTalliesResultsTypes<T+1>();
    }

    // ++ Global results maps
      
    template <typename... TallyTypes>
    struct GlobalResultsMapTupleBuilder{
      using type = void;
    };
    
    template <typename... TallyTypes>
    struct GlobalResultsMapTupleBuilder<std::tuple<TallyTypes...>>{
      using type = std::tuple<typename TallyTypes::ResultsMapType...>;
    };

    template <typename... TallyTypes>
    using GlobalResultsMapType = typename GlobalResultsMapTupleBuilder<TallyTypes...>::type;
    
    
    class Results{
      
    public:
      using Maps = GlobalResultsMapType<typesGenericTallies>;

    private:
      Maps maps;

      // ++ Extract results functions
      template<size_t I>
      typename std::enable_if<I >= std::tuple_size<typesGenericTallies>::value, void>::type      
      extractResults(pen_commonTallyCluster&, const unsigned, const unsigned long long){}
      
      template<size_t I = 0>
      typename std::enable_if<I < std::tuple_size<typesGenericTallies>::value, void>::type      
      extractResults(pen_commonTallyCluster& tallies, const unsigned iTally, const unsigned long long nhists){

	using TallyType = typename std::tuple_element<I, typesGenericTallies>::type;	
	
	if(tallies.isCreated<TallyType>(iTally)){
	  //The tally type match, get the name
	  std::string tallyName = tallies.getTallyName(iTally);

	  //Save the results
	  std::get<I>(maps)[tallyName] = tallies.getResults<TallyType>(iTally, nhists);
	}
	else{
	  //Try with the next tally type
	  extractResults<I+1>(tallies, iTally, nhists);
	}
      }

      // ++ Clear functions
      template<size_t I>
      typename std::enable_if<I >= std::tuple_size<typesGenericTallies>::value, void>::type      
      clearInner(){}
      
      template<size_t I = 0>
      typename std::enable_if<I < std::tuple_size<typesGenericTallies>::value, void>::type      
      clearInner(){
	std::get<I>(maps).clear();
	clearInner<I+1>();
      }
	
      // ++ Print functions

      // - Print results element
      template<class T>
      typename std::enable_if<measurements::is_results<T>::value, void>::type
      printResult(const T& v,
		  const char* filename,
		  const unsigned nSigma,
		  const bool coordinates,
		  const bool binNumber,
		  const bool onlyEffective) const {
	FILE* fout = fopen(filename, "w");
	if(fout != nullptr){
	  v.print(fout,nSigma,coordinates,binNumber,onlyEffective);
	  fclose(fout);
	}
      }
      
      template<class T>
      typename std::enable_if<std::is_arithmetic<T>::value, void>::type
      printResult(const std::vector<T>& v,
		  const char* filename,
		  const unsigned ,
		  const bool ,
		  const bool binNumber,
		  const bool ) const {
	FILE* fout = fopen(filename, "w");
	if(fout != nullptr){

	  if(std::is_floating_point<T>::value){

	    if(binNumber){
	      for(size_t i = 0; i < v.size(); ++i){
		fprintf(fout,"%4lu %15.5E\n",
			static_cast<unsigned long>(i),
			static_cast<double>(v[i]));
	      }
	    }
	    else{
	      for(size_t i = 0; i < v.size(); ++i){
		fprintf(fout,"%15.5E\n",
			static_cast<double>(v[i]));
	      }
	    }
	
	  }else{
	    if(std::is_signed<T>::value){
	  
	      if(binNumber){
		for(size_t i = 0; i < v.size(); ++i){
		  fprintf(fout,"%4lu %lld\n",
			  static_cast<unsigned long>(i),
			  static_cast<long long int>(v[i]));
		}
	      }
	      else{
		for(size_t i = 0; i < v.size(); ++i){
		  fprintf(fout,"%lld\n",
			  static_cast<long long int>(v[i]));
		}
	      }	      
	    }
	    else{

	      if(binNumber){
		for(size_t i = 0; i < v.size(); ++i){
		  fprintf(fout,"%4lu %lld\n",
			  static_cast<unsigned long>(i),
			  static_cast<long long unsigned>(v[i]));
		}
	      }
	      else{
		for(size_t i = 0; i < v.size(); ++i){
		  fprintf(fout,"%lld\n",
			  static_cast<long long unsigned>(v[i]));
		}
	      }	      
	    }
	  }	  
	  fclose(fout);
	}
      }

      template<class T>
      typename std::enable_if<!std::is_arithmetic<T>::value, void>::type
      printResult(const std::vector<T>& ,
		  const char* ,
		  const unsigned ,
		  const bool ,
		  const bool ,
		  const bool ) const {
	//Avoid printing non arithmetic types
      }

      // - Print all results for a specific tally type
      template<class TallyType, size_t I = 0>
      typename std::enable_if<I >= std::tuple_size<typename TallyType::ResultsTypes>::value, void>::type
      printTypeInner(const char*,
		     const unsigned ,
		     const bool ,
		     const bool ,
		     const bool ,
		     const typename TallyType::ResultsTypes&) const {}

      template<class TallyType, size_t I = 0>
      typename std::enable_if<I < std::tuple_size<typename TallyType::ResultsTypes>::value, void>::type
      printTypeInner(const char* prefix,
		     const unsigned nSigma,
		     const bool coordinates,
		     const bool binNumber,
		     const bool onlyEffective,
		     const typename TallyType::ResultsTypes& results) const {
	
	std::string filename(prefix);
	filename += "_";
	filename += TallyType::tallyID();
	filename += "_" + std::to_string(I) + ".dat";

	printResult(std::get<I>(results),
		    filename.c_str(),
		    nSigma,
		    coordinates,
		    binNumber,
		    onlyEffective);
	printTypeInner<TallyType, I+1>(prefix,
				       nSigma,
				       coordinates,
				       binNumber,
				       onlyEffective,
				       results);
      }

      // - Print all results for all tally types
      template<size_t I>
      typename std::enable_if<I >= std::tuple_size<typesGenericTallies>::value, void>::type      
      printInner(const char*,
		 const unsigned ,
		 const bool ,
		 const bool ,
		 const bool) const {}
      
      template<size_t I = 0>
      typename std::enable_if<I < std::tuple_size<typesGenericTallies>::value, void>::type      
      printInner(const char* prefix,
		 const unsigned nSigma,
		 const bool coordinates,
		 const bool binNumber,
		 const bool onlyEffective) const {

	using TallyType = typename std::tuple_element<I, typesGenericTallies>::type;	

	//Iterate over the results map for this tally type and write the results
	for(auto it = std::get<I>(maps).cbegin(); it != std::get<I>(maps).end(); ++it){
	  std::string filename(prefix);
	  filename += it->first;
	  printTypeInner<TallyType>(filename.c_str(),
				    nSigma,
				    coordinates,
				    binNumber,
				    onlyEffective,
				    it->second);
	}

	printInner<I+1>(prefix,
			nSigma,
			coordinates,
			binNumber,
			onlyEffective);
      }
      
    public:

      template<size_t T>
      constexpr const auto& read() const{	
	static_assert(T < std::tuple_size<typesGenericTallies>::value,
		      "Invalid tally index. Unable to get tally results maps");
	
	return std::get<T>(maps);
      }
      
      template<size_t T>
      constexpr auto& get(){	
	static_assert(T < std::tuple_size<typesGenericTallies>::value,
		      "Invalid tally index. Unable to get tally results maps");
	
	return std::get<T>(maps);
      }
      
      template<class TallyType>
      constexpr auto getType(){
	
	//Get the tally index within the tuple
	constexpr const unsigned index = typeIndex<TallyType>();
	static_assert(index < std::tuple_size<typesGenericTallies>::value,
		      "Invalid tally type. Unable to get tally results maps");
	
	return std::get<index>(maps);
      }

      inline void update(pen_commonTallyCluster& tallies, const unsigned long long nhists){
	const unsigned nTallies = tallies.numTallies();

	//Extract all tallies' results
	for(unsigned iTally = 0; iTally < nTallies; ++iTally){
	  extractResults(tallies, iTally, nhists);
	}
      }

      inline bool update(pen_commonTallyCluster& tallies, const size_t iTally, const unsigned long long nhists){
	const unsigned nTallies = tallies.numTallies();

	if(iTally < nTallies){
	  extractResults(tallies, iTally, nhists);
	  return true;
	}
	return false;
      }
      
      inline bool update(pen_commonTallyCluster& tallies, const std::string tallyName, const unsigned long long nhists){
	int iTally = tallies.getTallyIndex(tallyName);
	if(iTally >= 0){
	  extractResults(tallies, iTally, nhists);
	  return true;
	}
	return false;
      }
      
      inline void print(const char* prefix,
			const unsigned nSigma,
			const bool coordinates,
			const bool binNumber,
			const bool onlyEffective = false) const {
	printInner(prefix,nSigma,coordinates,binNumber,onlyEffective);
      }

      inline void clear(){
	clearInner();
      }
      
    };
  } //namespace tally
} // namespace penred

#endif
