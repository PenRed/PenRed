
//
//
//    Copyright (C) 2024-2025 Universitat de València - UV
//    Copyright (C) 2024-2025 Universitat Politècnica de València - UPV
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
//        vicent.gimenez.alventosa@gmail.com  (Vicent Giménez Alventosa)
//        sanolgi@upvnet.upv.es (Sandra Oliver Gil)
//    
//


#ifndef __PEN_SINGLES_TALLY__
#define __PEN_SINGLES_TALLY__

#include "pen_constants.hh"
#include <algorithm>
#include <cstdint>

class pen_Singles : public pen_genericTally<pen_particleState> {
  DECLARE_TALLY(pen_Singles,pen_particleState,SINGLES)

public:

  struct single{

    enum boolMaskEnum : uint8_t{
      SCATTERED            = 1 << 0, //The particle has been scattered outside detectors after production
      PRODUCED_IN_DETECTOR = 1 << 1, //The particle has been produced in a detector
      FIRST_GENERATION     = 1 << 2, //The particle belongs to the first genertion  (ILB[0]=1) 
      SECOND_GENERATION    = 1 << 3, //The particle belongs to the second genertion (ILB[0]=2)
      FROM_ANNIHILATION    = 1 << 4, //The particle has been created from a positron annihilation
      PILEUP               = 1 << 5, //The pulse mixes hits from different histories due pileup
    };

    static constexpr const size_t maxBuffSize = 200;
    static constexpr const size_t dataSize =
      sizeof(float)*5 + sizeof(double) + sizeof(uint8_t)*3 + sizeof(unsigned long long);
    
    double E, x, y, z, weight, t;
    std::array<uint8_t, 3> info; // (mask, kpar, interaction)
    unsigned long long hist;

    constexpr single() noexcept : E(0.0), x(0.0), y(0.0), z(0.0),
      weight(0.0), t(0.0), info{0,0,0}, hist(0)
    {}
    constexpr single(const double de,
		     const double xIn, const double yIn, const double zIn,
		     const double tIn, const double w, const uint8_t boolMask,
		     const uint8_t kpar, const uint8_t interaction,
		     const unsigned long long histIn) noexcept :
      E(de),
      x(de*xIn),
      y(de*yIn),
      z(de*zIn),
      weight(de*w),
      t(tIn),
      info{boolMask, kpar, interaction},
      hist(histIn)
    {}

    //Bool mask functions
    static constexpr bool isScattered(const uint8_t mask) noexcept {return mask & SCATTERED;}
    static constexpr bool isProducedInDetector(const uint8_t mask) noexcept {return mask & PRODUCED_IN_DETECTOR;}
    static constexpr bool isFirstGeneration(const uint8_t mask) noexcept {return mask & FIRST_GENERATION;}
    static constexpr bool isSecondGeneration(const uint8_t mask) noexcept {return mask & SECOND_GENERATION;}
    static constexpr bool isFromAnnihilation(const uint8_t mask) noexcept {return mask & FROM_ANNIHILATION;}
    static constexpr bool isPileup(const uint8_t mask) noexcept {return mask & PILEUP;}

    constexpr bool isScattered() const noexcept {return isScattered(info[0]);}
    constexpr bool isProducedInDetector() const noexcept {return isProducedInDetector(info[0]);}
    constexpr bool isFirstGeneration() const noexcept {return isFirstGeneration(info[0]);}
    constexpr bool isSecondGeneration() const noexcept {return isSecondGeneration(info[0]);}
    constexpr bool isFromAnnihilation() const noexcept {return isFromAnnihilation(info[0]);}
    constexpr bool isPileup() const noexcept {return isPileup(info[0]);}

    inline std::string stringify() const noexcept {
      char auxBuff[maxBuffSize];
      toBufferFinal(auxBuff, maxBuffSize);
      std::string ret(auxBuff);
      ret.pop_back();
      return ret;
    }

    inline void add(const single& s) noexcept {

      E += s.E;
      x += s.x;
      y += s.y;
      z += s.z;
      weight += s.weight;
    }

    inline void addEnergy(const double de, const double w) noexcept {
      //Scale energy increment according to weight ratio
      const double origW = weight/E;
      const vector3D<double> porig = pos();
      E = E + de*w/origW;

      //Correct parameters
      x = porig.x*E;
      y = porig.y*E;
      z = porig.z*E;
      weight = origW*E;
      
    }

    inline vector3D<double> pos() const noexcept {
      return vector3D<double>(x/E,y/E,z/E);
    }

    inline double dist(const single& s) const noexcept{
      return pos().dist(s.pos());
    }

    inline int toBuffer(char* b, size_t max) const {
      return snprintf(b, max, "%15.5E %15.5E %15.5E %15.5E "
		      "%25.15E %15.5E %u %u %u %llu\n",
		      E, x, y, z, weight, t,
		      info[0], info[1], info[2], hist);
    }

    inline int toBufferFinal(char* b, size_t max) const {
      return snprintf(b, max, "%15.5E %15.5E %15.5E %15.5E "
		      "%15.5E %25.15E %u %u %u %llu\n",
		      E, x/E, y/E, z/E, weight/E, t,
		      info[0], info[1], info[2], hist);
    }

    inline void toBufferB(unsigned char* b, size_t& pos){
      float auxf[] = {(float)E, (float)x, (float)y, (float)z, (float)weight};
      memcpy(&b[pos], auxf, 5*sizeof(float));
      pos += sizeof(float)*5;

      memcpy(&b[pos], &t, sizeof(double));
      pos += sizeof(double);      
      memcpy(&b[pos], info.data(), 3*sizeof(uint8_t));
      pos += 3*sizeof(uint8_t);
      memcpy(&b[pos], &hist, sizeof(unsigned long long));
      pos += sizeof(unsigned long long);
    }

    inline void toBufferFinalB(unsigned char* b, size_t& pos){
      float auxf[] = {(float)E, float(x/E), float(y/E), float(z/E), float(weight/E)};
      memcpy(&b[pos], auxf, 5*sizeof(float));
      pos += sizeof(float)*5;

      memcpy(&b[pos], &t, sizeof(double));
      pos += sizeof(double);
      memcpy(&b[pos], info.data(), 3*sizeof(uint8_t));
      pos += 3*sizeof(uint8_t);
      memcpy(&b[pos], &hist, sizeof(unsigned long long));
      pos += sizeof(unsigned long long);
    }
    
    inline bool read(FILE* f, unsigned long& offset){

      //Set the file position
      std::fseek(f, offset, SEEK_SET);

      offset += dataSize;

      float auxf[5];
      fread(static_cast<void*>(auxf), sizeof(float), 5, f);
      E = auxf[0];
      x = auxf[1];
      y = auxf[2];
      z = auxf[3];
      weight = auxf[4];

      fread(static_cast<void*>(&t), sizeof(double), 1, f);
      fread(static_cast<void*>(info.data()), sizeof(uint8_t), 3, f);
      if(fread(static_cast<void*>(&hist), sizeof(unsigned long long), 1, f) == 1)
	return true;
      else
	return false;
      
    }
    
    inline bool operator<(const single& o) const noexcept{
      if(t == o.t){
	return E > o.E; //Set first high or positive energy
      }
      return t < o.t;
    }
  };

  struct singlesBuffer{
  public:
    static constexpr const size_t baseSize = 200000;
  private:
    std::vector<single> buffer;
    size_t n;
    size_t nflushes;
  public:
    singlesBuffer() noexcept : buffer(baseSize), n(0), nflushes(0) {}
    inline void store(const double de,
		      const double x, const double y, const double z,
		      const double t, const double w, const uint8_t boolMask,
		      const uint8_t kpar, const uint8_t interaction,
		      const unsigned long long hist) noexcept {
      
      if(n >= buffer.size()){
	buffer.resize(buffer.size() + baseSize/10);
      }
      buffer[n++] = single(de, x, y, z, t, w, boolMask, kpar, interaction, hist);
    }

    void reduce(const size_t start, const bool saveScatter,
		const double dt, const double tmin, const double tmax);
    
    inline std::vector<unsigned char> flush(){

      if(n == 0)
	return std::vector<unsigned char>();

      //Sort data
      std::sort(buffer.begin(), buffer.begin()+n);

      //Create a results buffer
      std::vector<unsigned char> result;
      result.resize(single::dataSize*n);

      //Save data in the results buffer
      size_t pos = 0;
      for(size_t i = 0; i < n; ++i){
	buffer[i].toBufferB(result.data(), pos);
      }

      //Reset buffer
      n = 0;

      //Increase number of flushes
      ++nflushes;
      
      //Return data
      return result;
    }
    inline size_t flushes() const { return nflushes; }
    inline size_t size() const { return n; }
  };

private:
  
  double singleEmin, singleEmax;
  double tmin, tmax;
  double dt;

  double joinTime;

  bool pileup;
  bool scatter;

  std::array<unsigned long, constants::nParTypes> nInStack;

  //The tally will use one buffer per detector
  std::vector<singlesBuffer> buffers;
  //Save last history start position in each buffer
  std::vector<size_t> lastHistStart;
  //Save enabled detector indexes
  std::vector<unsigned> detInternalIndex;
  //Save sensible detectors
  std::vector<bool> kdets;
  //Save paths to information files
  std::vector<std::string> fInfoPaths;
  //Save paths to data files
  std::vector<std::string> fDataPaths;
  //Save data files offsets (Bytes from file beggining)
  std::vector<unsigned long> offsets;

  uint8_t actualMask;
  unsigned actualKdet;
  bool toDetect;
  bool simFinished;
  bool removeOnEnd;
  bool binary;
  bool skipBeginPart;
  unsigned long long lastHist;

  // Store last interaction index. Last values are for special situations:
  // 255 Particle just created and absorbed due local body eabs
  // 254 Particle absorbed just after entering a new material or body due local eabs
  uint8_t lastICol;

  const wrapper_geometry* geo;
  
public:

  pen_Singles() : pen_genericTally( USE_LOCALEDEP |
				    USE_BEGINPART |
				    USE_SAMPLEDPART |
				    USE_STEP |
				    USE_MOVE2GEO |
				    USE_INTERFCROSS |
				    USE_KNOCK |
				    USE_ENDHIST |
				    USE_ENDSIM)
  {}

  inline bool activeDet(){
    //Check detector
    if(actualKdet >= kdets.size())
      return false;
    else
      return kdets[actualKdet];
  }

  inline void count(const double de,
		    const double x, const double y, const double z,
		    const double t, const double w, const pen_KPAR kpar,
		    const unsigned long long nhist){
    
    if(!toDetect){
      if(isRealInteraction(kpar, lastICol))
	actualMask |= single::SCATTERED;
      return;
    }
    
    if(de == 0.0)
      return;

    //Get the buffer
    singlesBuffer& buffer = buffers[detInternalIndex[actualKdet]];

    //Store single
    buffer.store(de, x, y, z, t, w, actualMask, static_cast<uint8_t>(kpar), lastICol, nhist);
  }
  
  void tally_localEdep(const unsigned long long nhist,
		       const pen_KPAR /*kpar*/,
		       const pen_particleState& state,
		       const double dE);
  
  void tally_beginPart(const unsigned long long nhist,
		       const unsigned kdet,
		       const pen_KPAR /*kpar*/,
		       const pen_particleState& state);

  void tally_sampledPart(const unsigned long long nhist,
			 const unsigned long long /*dhist*/,
			 const unsigned kdet,
			 const pen_KPAR /*kpar*/,
			 const pen_particleState& state);
  
  void tally_step(const unsigned long long nhist,
		  const pen_KPAR /*kpar*/,
		  const pen_particleState& state,
		  const tally_StepData& stepData);

  void tally_move2geo(const unsigned long long nhist,
		      const unsigned kdet,
		      const pen_KPAR /*kpar*/,
		      const pen_particleState& state,
		      const double /*dsef*/,
		      const double /*dstot*/);  

  void tally_interfCross(const unsigned long long /*nhist*/,
			 const unsigned kdet,
			 const pen_KPAR /*kpar*/,
			 const pen_particleState& /*state*/);

  inline void tally_knock(const unsigned long long /*nhist*/,
			  const pen_KPAR /*kpar*/,
			  const pen_particleState& /*state*/,
			  const int icol){
    //Update last interaction
    lastICol = static_cast<uint8_t>(icol);
  }

  void tally_endHist(const unsigned long long /*nhist*/);

  inline void tally_endSim(const unsigned long long /*nhist*/){
    simFinished = true;
  }

  int configure(const wrapper_geometry& geometry,
		const abc_material* const /*materials*/[pen_geoconst::NB],     
		const pen_parserSection& /*config*/, const unsigned verbose);
  void flush(const unsigned det);
  void flush();
  void saveData(const unsigned long long nhist) const;
  int sumTally(const pen_Singles& tally);

  void orderDetectorData(const unsigned idet, FILE* fout) const;
};

//Tally configuration reader
class tallyReader_Singles : public pen_configReader<tallyReader_Singles>{

public:
  
  enum errors{
    SUCCESS = 0,
    UNHANDLED = 1,
  };

  bool pileup;
  bool scatter;
  
  bool removeOnEnd;
  bool binary;
  
  double singleEmin, singleEmax;
  double tmin, tmax, tjoin;

  std::vector<unsigned> kdets;

  int storeElement(const std::string& pathInSection,
		   const pen_parserData& element,
		   const unsigned verbose);
  
  int beginArray(const std::string& pathInSection,
		 const size_t size,
		 const unsigned verbose);

  int storeArrayElement(const std::string& pathInSection,
			const pen_parserData& element,
			const size_t pos,
			const unsigned verbose);

  inline int endArray(const unsigned /*verbose*/) { return errors::SUCCESS; }
  
};

template<>
struct pen_format<tallyReader_Singles>{
  //By default, no format is defined
  static constexpr const char* format = R"===(

# Tally "Singles" reader configuration

reader-description "Tally to register singles reaching a set of detectors"

# Detectors
detectors/reader-description "Sensible detectors"
detectors/reader-value [1]
detectors/reader-required/type "required"

# Thresholds

## Energy
energy/single/min/reader-description "Minimum single energy to be stored, in eV"
energy/single/min/reader-value 0.0
energy/single/min/reader-required/type "required_if_exist"
energy/single/min/reader-required/value "energy/single/max"
energy/single/min/reader-conditions/lt/type "lesser"
energy/single/min/reader-conditions/lt/value "energy/single/max"
energy/single/min/reader-conditions/positive/type "positive"

energy/single/max/reader-description "Maximum single energy to be stored, in eV"
energy/single/max/reader-value 1.0e30
energy/single/max/reader-required/type "required_if_exist"
energy/single/max/reader-required/value "energy/single/min"
energy/single/max/reader-conditions/gt/type "greater"
energy/single/max/reader-conditions/gt/value "energy/single/min"

## Time
time/min/reader-description "Minimum time to be tallied in seconds"
time/min/reader-value 0.0
time/min/reader-required/type "required_if_exist"
time/min/reader-required/value "time/max"
time/min/reader-conditions/lt/type "lesser"
time/min/reader-conditions/lt/value "time/max"
time/min/reader-conditions/positive/type "positive"

time/max/reader-description "Maximum time to be tallied in seconds"
time/max/reader-value 1.0e30
time/max/reader-required/type "required_if_exist"
time/max/reader-required/value "time/min"
time/max/reader-conditions/gt/type "greater"
time/max/reader-conditions/gt/value "time/min"

time/join/reader-description "Minimum required time between singles signals, in seconds. Singals in the same time window will be added."
time/join/reader-value 1.0e-9
time/join/reader-required/type "optional"
time/join/reader-conditions/positive/type "positive"

## Events control
pileup/reader-description "Enable/disable pileup when joining events. If disabled, only hits proceding from the same history can be joined in a single."
pileup/reader-value true
pileup/reader-required/type "optional"

scatter/reader-description "Enable/disable saving scattered events."
scatter/reader-value true
scatter/reader-required/type "optional"

## File control
clear/reader-description "Enable/disable removing auxiliary data files after final data processing. Notice that those files are required to resume a simulation from a dump file."
clear/reader-value true
clear/reader-required/type "optional"

binary/reader-description "Enable/disable binary output."
binary/reader-value true
binary/reader-required/type "optional"
)===";
};


#endif
