
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

class pen_Singles : public pen_genericTally<pen_particleState> {
  DECLARE_TALLY(pen_Singles,pen_particleState)

private:

  struct single{

    static constexpr const size_t maxBuffSize = 150;
    static constexpr const size_t dataSize = sizeof(float)*5 + sizeof(double) + sizeof(unsigned long long);
    
    double E, x, y, z, t, weight;
    unsigned long long hist;

    constexpr single() : E(0.0), x(0.0), y(0.0), z(0.0), t(0.0),
			 weight(0.0), hist(0)
    {}
    inline single(const double de,
		  const double xIn, const double yIn, const double zIn,
		  const double tIn, const double w, const unsigned long long histIn) noexcept :
      E(de),
      x(de*xIn),
      y(de*yIn),
      z(de*zIn),
      t(tIn),
      weight(de*w),
      hist(histIn)
    {}

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
		      "%25.15E %15.5E %llu\n",
		      E, x, y, z, t, weight, hist);
    }

    inline int toBufferFinal(char* b, size_t max) const {
      return snprintf(b, max, "%15.5E %15.5E %15.5E %15.5E "
		      "%25.15E %15.5E %llu\n",
		      E, x/E, y/E, z/E, t, weight/E, hist);
    }

    inline void toBufferB(unsigned char* b, size_t& pos){
      float auxf[] = {(float)E, (float)x, (float)y, (float)z, (float)weight};
      memcpy(&b[pos], auxf, 5*sizeof(float));
      pos += sizeof(float)*5;

      memcpy(&b[pos], &t, sizeof(double));
      pos += sizeof(double);      
      memcpy(&b[pos], &hist, sizeof(unsigned long long));
      pos += sizeof(unsigned long long);
    }

    inline void toBufferFinalB(unsigned char* b, size_t& pos){
      float auxf[] = {(float)E, float(x/E), float(y/E), float(z/E), float(weight/E)};
      memcpy(&b[pos], auxf, 5*sizeof(float));
      pos += sizeof(float)*5;

      memcpy(&b[pos], &t, sizeof(double));
      pos += sizeof(double);      
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
    static constexpr const size_t baseSize = 200000;
  private:
    std::vector<single> buffer;
    size_t n;
    size_t nflushes;
  public:
    singlesBuffer() noexcept : buffer(baseSize), n(0), nflushes(0) {}
    inline void store(const double de,
		      const double x, const double y, const double z,
		      const double t, const double w,
		      const unsigned long long hist) noexcept {
      
      if(n < buffer.size()){
	buffer[n++] = single(de, x, y, z, t, w, hist);
      }else{
	buffer.resize(buffer.size() + baseSize/10);
	buffer[n++] = single(de, x, y, z, t, w, hist);
      }
    }

    inline void reduce(const size_t start, const double dt, const double tmin, const double tmax){

      if(n <= start+1)
	return;      

      //std::vector<single> auxBuff(buffer.cbegin() + start, buffer.cbegin() + n);
      
      //Sort specified data range
      std::sort(buffer.begin() + start, buffer.begin() + n);

      //std::vector<single> auxBuffSort(buffer.cbegin() + start, buffer.cbegin() + n);

      // Add singles within the same time window
      
      //Create auxiliary single
      single auxSing;
      size_t iAux = start;
      
      //Get the first single with positive energy
      size_t ifirst = iAux;
      auxSing = buffer[ifirst];
      while(auxSing.E < 0.0 && ifirst < n-1){
	//Get the next one
	auxSing = buffer[++ifirst]; 
      }
      
      if(ifirst != iAux){
	//The first single with positive energy is not at the first position, correct it
	buffer[ifirst] = buffer[iAux];
      }

      if(auxSing.E < 0.0){
	penred::logs::logger::printf(penred::logs::SIMULATION,
				     "No non negative single found!\n"
				     "   First Energy: %15.5E",
				     buffer[ifirst].E);
      }

      //double lastE = auxSing.E;
      for(size_t i = iAux+1; i < n; ++i){
	if(buffer[i].E < 0.0){ //Check if the next pulse has a negative energy 	  
	  auxSing.addEnergy(buffer[i].E, buffer[i].weight/buffer[i].E);
	}
	else if(buffer[i].t - auxSing.t < dt){ //Check time window
	  auxSing.add(buffer[i]);	  
	}else{
	  //Pulse with positive energy outside the join window

	  //lastE = auxSing.E;

	  //Save the current single
	  if(auxSing.E > 0.0){
	    if(auxSing.t >= tmin && auxSing.t <= tmax)
	      buffer[iAux++] = auxSing;
	  }
	  /*
	  else if(auxSing.E <= -1.0e-6)
	    penred::logs::logger::printf(penred::logs::SIMULATION,
					 "Warning: Single with negative energy (%E eV)."
					 " Previous single energy: %E eV\n"
					 "       i first: %lu\n"					 
					 "         Start: %lu\n"
					 "          iAux: %lu\n"
					 "             n: %lu\n",
					 auxSing.E, lastE, ifirst, start, iAux, n);
	  */
	  
	  //Get the next one
	  auxSing = buffer[i];
	}
	
      }
      
      //Save the last single
      if(auxSing.E > 0.0){
	buffer[iAux++] = auxSing;
      }
      /*
      else if(auxSing.E <= -1.0e-6){

	std::string text;
	for(size_t k = 0; k < auxBuff.size(); ++k){
	  text.append(auxBuffSort[k].stringify() + "  |  " + auxBuff[k].stringify());
	  text.append("\n");
	}
	
	penred::logs::logger::printf(penred::logs::SIMULATION,
				     "Warning: Final single with negative energy (%E eV).\n"
				     "         Previous single energy: %E eV\n"
				     " Previous history single energy: %E eV\n"
				     "       i first: %lu\n"					 
				     "         Start: %lu\n"
				     "          iAux: %lu\n"
				     "             n: %lu\n"
				     "         tinit: %.25E\n"				     
				     "         tlast: %.25E\n"
				     "         Einit: %.15E\n"				     
				     "         Elast: %.15E\n"
				     " Singles list (%lu):\n%s\n",
				     auxSing.E, lastE, buffer[start > 0 ? start-1 : 0].E, ifirst, start, iAux, n,
				     auxSing.t, buffer[n-1].t, buffer[iAux].E, buffer[n-1].E,
				     auxBuff.size(), text.c_str());
      }
      */

      
      //Update the number of elements in the buffer
      n = iAux;      
    }
    
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
  
  double singleEmin, singleEmax;
  double tmin, tmax;
  double dt;

  double joinTime;

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

  unsigned actualKdet;
  bool toDetect;
  bool simFinished;
  bool removeOnEnd;
  bool binary;
  bool skipBeginPart;

  const wrapper_geometry* geo;
  
public:

  pen_Singles() : pen_genericTally( USE_LOCALEDEP |
				    USE_BEGINPART |
				    USE_SAMPLEDPART |
				    USE_STEP |
				    USE_MOVE2GEO |
				    USE_INTERFCROSS |
				    USE_ENDHIST)
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
		    const double t, const double w,
		    const unsigned long long nhist){

    if(!toDetect)
      return;
    
    if(de == 0.0)
      return;

    //Get the buffer
    singlesBuffer& buffer = buffers[detInternalIndex[actualKdet]];
    
    //Store single
    buffer.store(de, x, y, z, t, w, nhist);
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
