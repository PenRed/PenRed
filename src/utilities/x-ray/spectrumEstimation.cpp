//
//
//    Copyright (C) 2023-2024 Universitat de València - UV
//    Copyright (C) 2023-2024 Universitat Politècnica de València - UPV
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

#include <random>

#include "PenRed.hh"
#include "pen_samplers.hh"
#include "pen_tallies.hh"
#include "pen_geometries.hh"

template <class T>
struct measurment{

protected:
  std::vector<T> aux;
  std::vector<T> v;
  std::vector<T> v2;
  unsigned long long nCP;
public:

  measurment() = default;
  
  measurment(const size_t l){
    resize(l);
  }
  
  virtual inline void resize(const size_t l){
    aux.resize(l);
    v.resize(l);
    v2.resize(l);

    reset();
  }

  inline void add(const size_t i, const T val){
    if(i < v.size())
      aux[i] += val;
  }

  inline int sum(const measurment<T>& m){

    if(v.size() != m.v.size())
      return -1;
    
    for(size_t i = 0; i < v.size(); ++i){
      v[i] += m.v[i];
    }

    for(size_t i = 0; i < v.size(); ++i){
      v2[i] += m.v2[i];
    }

    nCP += m.nCP;

    return 0;
  }

  inline int set(std::vector<T> vIn, std::vector<T> v2In, const unsigned long long nCPIn){
    if(vIn.size() != v2In.size() || v.size() != vIn.size())
      return -1;

    v.swap(vIn);
    v2.swap(v2In);
    nCP = nCPIn;
    
    return 0;
  }

  inline void checkPoint(){

    for(size_t i = 0; i < v.size(); ++i){
      v[i] += aux[i];
      v2[i] += aux[i]*aux[i];
    }

    std::fill(aux.begin(), aux.end(), static_cast<T>(0.0));

    nCP += 1;
  }

  inline std::pair<std::vector<T>, std::vector<T>> results() const{

    std::pair<std::vector<T>, std::vector<T>> res(v.size(),v.size());

    double iNCP = 1.0/static_cast<double>(nCP);
    for(size_t i = 0; i < v.size(); ++i){
      double q = v[i]*iNCP;
      double q2 = v2[i]*iNCP;
      double sig = (q2 - q*q)*iNCP;
      if(sig > 0.0)
	sig = sqrt(sig);
      else
	sig = 0.0;

      res.first[i] = q;
      res.second[i] = sig;
    }

    return res;
  }

  inline std::vector<T> sigma() const{

    std::vector<T> sigmas(v.size());
    std::fill(sigmas.begin(),
	      sigmas.end(),
	      static_cast<T>(0.0));

    double iNCP = 1.0/static_cast<double>(nCP);
    for(size_t i = 0; i < v.size(); ++i){
      double q = v[i]*iNCP;
      double q2 = v2[i]*iNCP;
      double sig = (q2 - q*q)*iNCP;
      if(sig > 0.0)
	sig = sqrt(sig);
      else
	sig = 0.0;

      sigmas[i] = sig;
    }

    return sigmas;
  }
  
  inline double maxAbsErel(const double absThreshold = 1.0) const{

    //Calculate threshold
    double threshold = -1.0;
    if(absThreshold < 1.0){
      double maxVal = *std::max_element(v.begin(), v.end(),
					[](const double a, const double b){
					  return std::fabs(a) < std::fabs(b);
					});
      threshold = maxVal*absThreshold;
    }
    
    double max = 0.0;
    double iNCP = 1.0/static_cast<double>(nCP);
    for(size_t i = 0; i < v.size(); ++i){
      if(std::fabs(v[i]) > threshold){
	double q = v[i]*iNCP;
	double q2 = v2[i]*iNCP;
	double sig = (q2 - q*q)*iNCP;
	if(sig > 0.0){
	  sig = sqrt(sig);
	  double eRel = std::fabs(sig/q);
	  if(eRel > max)
	    max = eRel;
	}
      }
    }
    return max;
  }

  inline double meanErel() const{

    double mean = 0.0;
    unsigned long nonZero = 0;
    double iNCP = 1.0/static_cast<double>(nCP);
    for(size_t i = 0; i < v.size(); ++i){
      double q = v[i]*iNCP;
      double q2 = v2[i]*iNCP;
      double sig = (q2 - q*q)*iNCP;
      if(sig > 0.0){
	sig = sqrt(sig);
	mean += sig/q;
	++nonZero;
      }
    }
    return mean/static_cast<double>(nonZero);
  }

  virtual inline void clear(){
    aux.clear();
    v.clear();
    v2.clear();

    nCP = 0;
  }

  virtual inline void reset(){
    std::fill(aux.begin(),
	      aux.end(),
	      static_cast<T>(0.0));
    std::fill(v.begin(),
	      v.end(),
	      static_cast<T>(0.0));
    std::fill(v2.begin(),
	      v2.end(),
	      static_cast<T>(0.0));

    nCP = 0;
  }
};

template <class T>
struct measurmentRegular : measurment<T>{

  T dx;
  T minx;

  inline void add(T x, const T y){
    x -= minx;
    if(!std::signbit(x)){
      size_t i = x/dx;
      measurment<T>::add(i, y);
    }
  }
};

template <class T>
struct measurmentIrregular : measurment<T>{

  std::vector<T> xLow;
  std::vector<T> dx;

public:

  inline void add(const T x, const T y){

    if(x < xLow[0]){
      return;
    }
    
    for(size_t i = 1; i < xLow.size(); ++i){
      if(x < xLow[i]){
	measurment<T>::add(i-1, y);
	return;
      }
    }

    if(x < xLow.back() + dx.back()){
      measurment<T>::add(xLow.size()-1, y);
    }
  }

  inline void resize(const size_t l){
    xLow.resize(l);
    dx.resize(l);

    measurment<T>::resize(l);
  }

  inline void clear(){
    xLow.clear();
    dx.clear();

    measurment<T>::clear();
  }

  
  
}; 

struct spectrumStruct{
  double emin, dE;
  std::vector<double> weight;
  std::vector<double> eWeight;
  std::vector<double> cummulative;

  inline int updateCummulative(){
    
    cummulative[0] = 0.0;
    for(size_t i = 1; i < cummulative.size(); ++i){
      cummulative[i] = cummulative[i-1] + weight[i-1];
    }

    //Normalize
    double sum = cummulative.back();

    //Check if the spectrum is empty
    if(sum == 0.0)
      return -1;
    
    for(double& c : cummulative){
      c /= sum;
    }

    return 0;
  }

  inline double sample(pen_rand& random) const {

    // Get random number
    const double rand = random.rand();

    // Get interval
    const unsigned interval = seeki(&cummulative.front(),rand,cummulative.size());

    // Calculate sampled energy
    const double ElowInterval = emin + dE*static_cast<double>(interval);
    return ElowInterval+(rand-cummulative[interval])*dE;
    
  }

  inline void clearWeight(){
    std::fill(weight.begin(), weight.end(), 0.0);
    std::fill(eWeight.begin(), eWeight.end(), 0.0);
    std::fill(cummulative.begin(), cummulative.end(), 0.0);
  }

  inline void resize(const size_t l,
		     const double EminIn,
		     const double dEIn){
    emin = EminIn;
    dE = dEIn;

    weight.resize(l);
    eWeight.resize(l);
    cummulative.resize(l+1);
  }

  inline void clear(){
    weight.clear();
    eWeight.clear();
    cummulative.clear();
  }

  inline void fprint(FILE* fout, const bool normalize = true) const {

    fprintf(fout, "#    Elow (eV)      prob (1/eV)      err (1/eV)     cummulative\n");
    double norm = dE;
    if(normalize){
      double sum = std::accumulate(weight.begin(), weight.end(), 0.0);
      norm *= sum;
    }

    if(norm > 0.0){
      for(size_t i = 0; i < weight.size(); ++i){
	double Elow = emin + dE*(static_cast<double>(i) + 0.5);
	fprintf(fout, "%15.5E %15.5E %15.5E %15.5E\n",
		Elow, weight[i]/norm, eWeight[i]/norm, cummulative[i+1]);
      }
    }
    
  }

  inline void print(const bool normalize = true) const {

    double norm = dE;
    if(normalize){
      double sum = std::accumulate(weight.begin(), weight.end(), 0.0);
      norm *= sum;
    }

    for(size_t i = 0; i < weight.size(); ++i){
      double Elow = emin + dE*(static_cast<double>(i) + 0.5);
      printf("%15.5E %15.5E %15.5E %15.5E\n",
	     Elow, weight[i]/norm, eWeight[i]/norm, cummulative[i+1]);
    }    
  }

};

struct objectiveMu {
  std::string name;
  double mu;
  double meanZ;
  unsigned imat;  //Saves material index
  
  objectiveMu() = default;
  objectiveMu(const std::string nameIn,
	      const double muIn,
	      const double meanZIn,
	      const unsigned imatIn) : name(nameIn),
				       mu(muIn),
				       meanZ(meanZIn),
				       imat(imatIn){}
};

struct filter{
  double thickness;
  double meanZ;
  unsigned imat;  //Saves material index (0 is void)

  filter() = default;
  filter(const double thicknessIn,
	 const double meanZIn,
	 const unsigned imatIn) : thickness(thicknessIn),
				  meanZ(meanZIn),
				  imat(imatIn){}

};

struct systemGeo{

private:

  //Know values
  double detectorLength, detectorLength05;
  double detectorWidth, detectorWidth05;
  double distSource2Det;
  
  //Unknown values
  double rInCol2, rOutCol2; //Inlet and outlet collimator radious
public:
  
  //Unknown values
  double distPSF2Filter;
  double dist2Col;
  double hCol;

  systemGeo() : detectorLength(1.0e35),
		detectorLength05(1.0e35),
		detectorWidth(1.0e35),
		detectorWidth05(1.0e35),
		distSource2Det(50.0),
		rInCol2(1.2),
		rOutCol2(4.0),
		distPSF2Filter(0.0),
		dist2Col(1.0),
		hCol(4.5){}

  //Set functions
  inline void setDetectorLength(const double v){
    detectorLength = v;
    detectorLength05 = v/2.0;
  }

  inline void setDetectorWidth(const double v){
    detectorWidth = v;
    detectorWidth05 = v/2.0;
  }

  inline void setDistSource2Det(const double v){
    distSource2Det = v;
  }
  
  inline void setRInCol(const double v){
    rInCol2 = v*v;
  }

  inline void setROutCol(const double v){
    rOutCol2 = v*v;
  }

  //Read functions
  inline double readRInCol() const { return sqrt(rInCol2);}
  inline double readRInCol2() const { return rInCol2;}

  inline double readROutCol() const { return sqrt(rOutCol2);}
  inline double readROutCol2() const { return rOutCol2;}

  inline double readDetectorLength() const { return detectorLength; }
  inline double readDetectorWidth() const { return detectorWidth; }
  inline double readDistSource2Det() const { return distSource2Det; }

  inline bool isDetected(const pen_particleState& state) const {

    //Calculate the particle position on collimator inlet and outlet 
    const double ds2Col = dist2Col/state.W;
    const double ds2Out = (dist2Col + hCol)/state.W;

    const double xOnCol = state.X + state.U*ds2Col;
    const double yOnCol = state.Y + state.V*ds2Col;

    const double xOnOut = state.X + state.U*ds2Out;
    const double yOnOut = state.Y + state.V*ds2Out;

    const double dOnCol2 = xOnCol*xOnCol + yOnCol*yOnCol;
    const double dOnOut2 = xOnOut*xOnOut + yOnOut*yOnOut;

    //Check if the particle is stopped by the collimator
    if(dOnCol2 > rInCol2 || dOnOut2 > rOutCol2){
      return false;
    }

    //The particle cross the collimator, check if it is detected
    const double ds2Det = (hCol + distSource2Det)/state.W;

    const double xOnDet = state.X + state.U*ds2Det;
    const double yOnDet = state.Y + state.V*ds2Det;

    if(xOnDet > detectorLength05 || xOnDet < -detectorLength05)
      return false;
    if(yOnDet > detectorWidth05 || yOnDet < -detectorWidth05)
      return false;

    return true;
  }
};

struct spectrumInspector{

private:
  std::vector<objectiveMu> objective;
  std::vector<size_t> addOrder;
  
  std::vector<double> muRhos;
  std::vector<double> eMuRhos;
  std::vector<double> diffMusRel;
  size_t maxDiff;
public:
  
  std::vector<filter> filters;
  systemGeo system;
  
  
  spectrumStruct spectrum;

  inline double relMaxDiff() const {return diffMusRel[maxDiff];}

  inline double foo() const {

    //Calculate the difference between the maximum discrepance and the other ones
    double maxRelDiffsBetweenDiffs = diffMusRel[maxDiff]/eMuRhos[0];
    size_t posOfMaxDiffsBetweenDiffs = 0;

    return 0.0;
    
  }
  
  inline size_t nFilters() const {return filters.size();}
  inline size_t nObjective() const {return muRhos.size();}

  inline void clearFilters(){
    filters.clear();
    std::fill(muRhos.begin(), muRhos.end(), 0.0);
    std::fill(eMuRhos.begin(), eMuRhos.end(), 1.0e35);    
    std::fill(diffMusRel.begin(), diffMusRel.end(), 1.0e35);
  }

  inline int setMu(const size_t imat, const double mu, const double sigma){

    bool found = false;
    size_t pos;
    for(size_t i = 0; i < objective.size(); ++i){
      if(objective[i].imat == imat){
	pos = i;
	found = true;
	break;
      }
    }

    if(!found)
      return -1;
    
    muRhos[pos] = mu;
    eMuRhos[pos] = 100.0*sigma/mu;
    double diff = mu/objective[pos].mu - 1.0;
    diffMusRel[pos] = diff;

    if(fabs(diff) > fabs(diffMusRel[maxDiff]))
      maxDiff = pos;

    return 0;
  }

  inline double sumDiffs2() const {
    return std::accumulate(diffMusRel.begin(),
			   diffMusRel.end(), 1.0,
			   [](const double a, const double b){
			     return a + b*b;
			   });
  }

  inline void addObjectiveMat(const unsigned imat,
			      const std::string name,
			      const double mu,
			      const double meanZ){

    //Create and add the new objective material
    objectiveMu toInsert(name, mu, meanZ, imat);
    toInsert.imat = imat;
    objective.push_back(toInsert);

    //Create results values for this material
    muRhos.push_back(0.0);
    eMuRhos.push_back(1.0e35);
    diffMusRel.push_back(1.0e35);

    //Reset maximum difference value
    maxDiff = 0;

    //Find the new objective position acording to Z value
    for(size_t i = 0; i < objective.size()-1; ++i){
      if(objective[i].meanZ > meanZ){

	//The new element must be inserted here
	
	//Move next elements
	for(size_t j = objective.size()-1; j > i; --j){
	  objective[j] = objective[j-1];
	}

	//Insert the new one
	objective[i] = toInsert;

	break;
      }
    }
  }

  inline void print(const bool printFilters = true,
		    const bool printMu = true,
		    const bool printSpectrum = true) const {

    if(printFilters){
      
      printf("\n    Filter materials: \n\n");
      printf(" Sim Mat Index      Z      Width(cm) \n");
      for(const filter& e : filters){
	printf("     %4u       %8.4f    %.4f\n",
	       e.imat, e.meanZ, e.thickness);
      }
    }

    if(printMu){
      printf("\n  Objective mu/rho values:\n\n");
      printf(" %-20s Index      Z      mu/rho(1/cm)     err(%%)     diff(%%)\n", "Name");
      for(unsigned i = 0; i < muRhos.size(); ++i){
	printf(" %-20s  %4u  %8.4f %13.5E      %.4f     %.4f\n",
	       objective[i].name.c_str(), objective[i].imat, objective[i].meanZ,
	       muRhos[i], eMuRhos[i], 100.0*diffMusRel[i]);
      }
    }

    if(printSpectrum){
      printf("\n\n");
      spectrum.print();
    }
  }

  inline void fprint(FILE* fout,
		     const bool printFilters = true,
		     const bool printMu = true,
		     const bool printSpectrum = true) const {

    if(printFilters){
      
      fprintf(fout,"\n#    Filter materials: \n\n");
      fprintf(fout,"# Sim Mat Index      Z      Width(cm) \n");
      for(const filter& e : filters){
	fprintf(fout,"#     %4u       %8.4f    %.4f\n",
		e.imat, e.meanZ, e.thickness);
      }
    }

    if(printMu){
      fprintf(fout,"\n# Objective mu/rho values:\n\n");
      fprintf(fout,"# %-20s index      Z      mu/rho(1/cm)     err(%%)     diff(%%)\n", "Name");
      for(unsigned i = 0; i < muRhos.size(); ++i){
	fprintf(fout,"# %-20s  %4u  %8.4f %13.5E      %.4f     %.4f\n",
		objective[i].name.c_str(), objective[i].imat, objective[i].meanZ,
		muRhos[i], eMuRhos[i], 100.0*diffMusRel[i]);
      }
    }

    if(printSpectrum){
      fprintf(fout, "\n\n");
      spectrum.fprint(fout);
    }
  }  

  inline void printRowHeader(FILE* fout) const {
    fprintf(fout,"#");
    for(const auto& e : filters){
      fprintf(fout," Filter mat %4u |",e.imat);
    }

    for(unsigned i = 0; i < muRhos.size(); ++i){
      if(i == muRhos.size()-1)
	fprintf(fout," Objective mat %4u",i);
      else
	fprintf(fout," Objective mat %4u %29s|",i," ");
    }

    fprintf(fout,"\n#");
    for(size_t i = 0; i < filters.size(); ++i)
      fprintf(fout,"    Width (cm)    ");
    for(size_t i = 0; i < muRhos.size(); ++i){
      fprintf(fout,"   Mu/rho(1/cm)     err(%%)     diff(%%)            ");
    }
    fprintf(fout,"\n");
  }
  
  inline void printRow(FILE* fout) const {

    for(size_t i = 0; i < filters.size(); ++i)
      fprintf(fout,"    %9.4f     ", filters[i].thickness);

    fprintf(fout," ");
    for(size_t i = 0; i < muRhos.size(); ++i){
      fprintf(fout,"%15.5E     %.4f    %+.2E %10s",
	      muRhos[i], eMuRhos[i], 100.0*diffMusRel[i]," ");
    }
    fprintf(fout,"\n");
  }
  
};

struct calcMuStruct{
  std::array<double, constants::MAXMAT> RANGET;
  std::array<double, constants::MAXMAT> RANGET2;

  unsigned long long nIter;

  inline void reset(const unsigned end){
    std::fill(RANGET.begin(),RANGET.begin() + end,0.0);
    std::fill(RANGET2.begin(),RANGET2.begin() + end,0.0);

    nIter = 0;
  }

  inline void sum(const calcMuStruct& b, const unsigned end){

    for(unsigned i = 0; i < end; ++i){
      RANGET[i] += b.RANGET[i];
    }

    for(unsigned i = 0; i < end; ++i){
      RANGET2[i] += b.RANGET2[i];
    }

    nIter += b.nIter;
  }
};


void simulate(const unsigned long long nIter,
	      const double tolerance,
	      const systemGeo* system,
	      const pen_context* pSimContext,
	      psfMemory_specificSampler* sampler,
	      const unsigned detMat,
	      measurmentRegular<double>* results,
	      int* seed1, int* seed2);

void calculateMu(const unsigned long long nIter,
		 const double tolerance,
		 const unsigned nMat,
		 const spectrumStruct* spectrum,
		 pen_context* contextObj,
		 calcMuStruct* results,
		 int* seed1, int* seed2);

template <class particleType>
bool move2geo(particleType& particle){

  pen_particleState& state = particle.getBaseState();
  //Check if particle has been generated at vacuum
  if(state.MAT == 0){
    //Move it to geometry system
    particle.jumpVolume();

    if(state.MAT == 0){
      return false;
    }
  }
  else{
    particle.updateMat();
    particle.updateBody();
  }
  return true;
}

template <class particleType>
inline bool absorb(particleType& particle,
		   pen_rand& penRand){

  const pen_particleState& state = particle.readBaseState();

  //Check if particle must be absorbed
  if(state.E < particle.getEABS()){ 

    //Check if the particle has already been annihilated
    if(state.E > 1.0e-6){
      // run annihilation process
      particle.annihilate(penRand);
    }
    return true;
  }

  //Particle can't be absorbed yet
  return false;
}

template <class particleType>
void simulatePart(particleType& particle,
		  pen_rand& randoms){

  pen_particleState& state = particle.getBaseState();

  //Check if particle scaped from geometry
  if(state.MAT == 0){
    return;
  }
  
  //Check if the particle must be absorbed
  if(absorb(particle,randoms))
    return;

  particle.START();
      
  for(;;){
    
    double ds, de;
    int icol;

    //Create step data structure
    tally_StepData stepData;
    
    //Calculate path to next iteration
    particle.JUMP(ds, randoms, particle.getDSMAX());
    
    //Move particle across geometry
    double lastEABS = particle.getEABS();
    particle.move(ds,stepData.softDE,
		  stepData.softX,stepData.softY,stepData.softZ,
		  randoms);    

    //Update step data
    stepData.update(particle);
    
    //Check if the particle must be absorbed at original material and body
    if(state.E < lastEABS){

      if(state.MAT > 0){
	// run annihilation process
	particle.annihilate(randoms);
      }
      //Set to the particle state the origin material and body
      state.IBODY = particle.lastBody();
      state.MAT = particle.lastMat();

      //Set position to the end of dsef
      double XL,YL,ZL;
      particle.lastPos(XL,YL,ZL);
      state.X = XL+state.U*particle.DSef();
      state.Y = YL+state.V*particle.DSef();
      state.Z = ZL+state.W*particle.DSef();
      return;
    }
      
    //Check if particle cross some interface
    if(particle.NCross() != 0){
      
      //Check if material has been changed
      if(particle.lastMat() != state.MAT){
	
	if(state.MAT == 0){
	  //New material is vacuum
	  return;
	}

	// Particle has been stopped at interface,
	// restart particle state
	particle.START();

	//Check if the particle must be absorbed
	if(absorb(particle,randoms))
	  return;	
      }
      else{
	// Particle has been stopped at interface,
	// restart particle state
	particle.START();

	//Check if the particle must be absorbed
	if(absorb(particle,randoms))
	  return;
      }
    }
    else{
      
      //No interface crossed, simulate interaction
      particle.KNOCK(de,icol,randoms);
      
      //Check if the particle must be absorbed
      if(absorb(particle,randoms)){
	return;
      }
    }
  }
}


int main(int argc, const char** argv){

  //Check the arguments
  if(argc < 2){
    printf("usage: %s configFile \n",argv[0]);
    return 1;
  }

  //Init random generator
  std::random_device rd;  // Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<double> rand01(0.0, 1.0);

  //**************************
  // Parse configuration
  //**************************

  //Parse configuration file
  pen_parserSection config;
  std::string errorLine;
  unsigned long errorLineNum;
  int err = parseFile(argv[1],config,errorLine,errorLineNum);
  
  if(err != INTDATA_SUCCESS){
    printf("Error parsing configuration.\n");
    printf("Error code: %d\n",err);
    printf("Error message: %s\n",pen_parserError(err));
    printf("Error located at line %lu, at text: %s\n",
	   errorLineNum,errorLine.c_str());
    return -1;
  }

  // Get verbose level
  //********************
  unsigned verbose = 2;
  int auxVerbose;
  if(config.read("simulation/verbose",auxVerbose) == INTDATA_SUCCESS){
    auxVerbose = std::max(0,auxVerbose);
    verbose = static_cast<unsigned>(auxVerbose);
    if(verbose > 1){
      printf("Verbose level set to %u\n",verbose);
    }
  }

  // Get minimum number of filters
  //********************************
  unsigned minFilters = 1;
  int auxMinFilters;
  if(config.read("simulation/minFilters",auxMinFilters) == INTDATA_SUCCESS){
    auxMinFilters = std::max(1,auxMinFilters);
    minFilters = static_cast<unsigned>(auxMinFilters);
    if(verbose > 1){
      printf("Minimum number of filters set to %u\n",minFilters);
    }
  }  
  
  // Histories per iteration
  //**************************
  double nHistsd;
  unsigned long long nHists;
  if(config.read("simulation/hists",nHistsd) == INTDATA_SUCCESS){
    nHists = static_cast<unsigned long long>(fabs(nHistsd));
  }else{
    nHists = 100000000;
  }

  // Solver iterations
  //**************************
  double nIterd;
  unsigned long long nIter;
  if(config.read("solver/iterations",nIterd) == INTDATA_SUCCESS){
    nIter = static_cast<unsigned long long>(fabs(nIterd));
  }else{
    nIter = 20;
  }
  
  // ** Number of threads
  //****************************
  unsigned nThreads = 1;
#ifdef _PEN_USE_THREADS_
  int nThreadsAux;
  err = config.read("simulation/threads", nThreadsAux);
  if(err != INTDATA_SUCCESS){
    nThreads = 1;
  }else{
    if(nThreadsAux <= 0){
      nThreads = std::max(static_cast<unsigned int>(2),
			  std::thread::hardware_concurrency());
    }
    else{
      nThreads = static_cast<unsigned>(nThreadsAux);
    }
  }
  
#endif

  // PSF sources
  //**************
  std::vector<psfMemory_specificSampler> psfs(nThreads);
  double Emax;
  
  //Read psf subsection
  pen_parserSection psfSection;
  err = config.readSubsection("psf", psfSection);
  if(err != INTDATA_SUCCESS){
    printf("Unable to read PSF section filename in "
	   "configuration file ('psf/')\n");
    return -3;
  }

  // Configure PSFs
  for(unsigned ith = 0; ith < nThreads; ++ith){

    psfMemory_specificSampler& psf = psfs[ith];

    psf.setThread(ith);
    double auxEmax;
    int errPsfConfig = psf.configure(auxEmax,
				     psfSection,
				     nThreads,
				     ith == 0 ? verbose : std::min(verbose,1u));
    if(errPsfConfig != 0){
      printf("Error configuring PSF source for thread %u.\n", ith);
      return -3;
    }

    if(ith == 0){
      Emax = auxEmax;
    }
  }

  //Once samplers have been configured for each thread, share thread 0 configuration
  for(unsigned ith = 1; ith < nThreads; ith++){
    psfMemory_specificSampler& psf = psfs[ith];
    int errSharing = psf.sharedConfig(psfs[0]);
    if(errSharing != 0){
      if(verbose > 0){
	printf("Error: Unable to share PSF source configuration "
	       "between thread 0 and %u\n",ith);
      }
      return -3;
    }
  }

  // ** Read parametrs for output spectrum
  //*********************************************

  //Read minimum output spectrum energy
  double outputEmin;
  err = config.read("spectrum/emin", outputEmin);
  if(err != INTDATA_SUCCESS){
    outputEmin = 1.0e3;
  }

  //Read maximum output spectrum energy
  double outputEmax;
  err = config.read("spectrum/emax", outputEmax);
  if(err != INTDATA_SUCCESS){
    outputEmax = Emax + 50;
  }

  //Read output spectrum number of bins
  int auxNbins;
  unsigned outputBins = 200;
  err = config.read("spectrum/nbins", auxNbins);
  if(err == INTDATA_SUCCESS && auxNbins > 0){
    outputBins = static_cast<unsigned>(auxNbins);
  }

  double outputDE = (outputEmax-outputEmin)/static_cast<double>(outputBins);
  
  if(verbose > 1){
    printf("Output spectrum energy parameters:\n");
    printf(" Minimum energy: %15.5E keV\n", outputEmin/1000.0);
    printf(" Maximum energy: %15.5E keV\n", outputEmax/1000.0);
    printf(" Number of bins: %5u\n", outputBins);
  }

  
  
  const double defEabs = outputEmin > 100.0 ? outputEmin-50.0 : 50.0;
  std::array<double, constants::nParTypes> minSimEabs;
  std::fill(minSimEabs.begin(), minSimEabs.end(), defEabs);

  // ** Initialization for each inspector and thread
  //**************************************************
  
  //Set initial seeds for each thread
  std::vector<int> seeds1(nThreads);
  std::vector<int> seeds2(nThreads);

  for(size_t ith = 0; ith < nThreads; ++ith){
    rand0(ith, seeds1[ith], seeds2[ith]);
  }

  //Init local results spectrums
  std::vector<measurmentRegular<double>> results(nThreads);
  //Get thread 0 results structure
  measurmentRegular<double>& globResults = results[0];

  //Create a vector for mu calculus results per thread
  std::vector<calcMuStruct> muThreads(nThreads);
  //Get thread 0 mu results structure
  calcMuStruct& globMuResults = muThreads[0];

  
  //Init results mesh for all threads
  for(size_t i = 0; i < results.size(); ++i){
    results[i].minx = outputEmin;
    results[i].dx = outputDE;
    results[i].resize(outputBins);
  }
  
  // ** Tolerance
  //****************
  double tolerance;
  err = config.read("simulation/tolerance", tolerance);
  if(err != INTDATA_SUCCESS){
    tolerance = 0.1;
  }else{
    tolerance = std::max(0.001,tolerance);
  }
  
  // ** Create simulation context  
  //*******************************

  //Create a context
  pen_context contextSim;

  // ** Materials 
  //***************
  std::string simMatFilePaths[constants::MAXMAT];  
  std::vector<std::string> filterMaterialsNames;
  err = config.ls("materials",filterMaterialsNames);
  if(err != INTDATA_SUCCESS){
    printf("Error: Section 'materials' not found in configuration file.\n");
    return -3;
  }

  if(filterMaterialsNames.size() < 1){
    printf("Error: No materials defined in 'materials' section.\n");
    return -3;
  }

  if(filterMaterialsNames.size() > constants::MAXMAT){
    printf("Error: Maximum number of materials exceeded.\n"
	   "   Maximum: %lu\n"
	   "   Used   : %lu\n",
	   static_cast<unsigned long>(constants::MAXMAT),
	   static_cast<unsigned long>(filterMaterialsNames.size()));
    return -3;
  }

  //Create materials in simulation context.
  //Create one extra dummy material for final detection
  int errmat = contextSim.setMats<pen_material>(filterMaterialsNames.size()+1);
  if(errmat != 0){
    printf("Error: Unable to create simulation context materials: %d.\n",errmat);
    return -3;
  }
  
  for(const std::string& matName : filterMaterialsNames){

    // Get filter path
    std::string matPath = std::string("materials/") + matName;
    pen_parserSection matSection;
    err = config.readSubsection(matPath, matSection);
    if(err != INTDATA_SUCCESS){
      printf("Error: Unable to read filter section '%s'.\n", matPath.c_str());
      return -3;
    }

    // ** Mat index
    int index;
    err = matSection.read("index",index);
    if(err != INTDATA_SUCCESS){
      printf("Error: Unable to read 'index' field in material '%s'\n",
	     matName.c_str());
      return -3;
    }

    if(index <= 0){
      printf("Error on material index assigned to material '%s' (%d). "
	     "Must be greater than zero\n",matName.c_str(), index);
      return -3;
    }

    unsigned imat = index-1;
    //Get material
    pen_material& mat = contextSim.getBaseMaterial(imat);
    
    //Set default C1, C2, WCC and WCR
    mat.C1=0.05;
    mat.C2=0.05;
    mat.WCC=defEabs;
    mat.WCR=defEabs;

    //Read absorption energies for each particle type
    for(unsigned ip = 0; ip < constants::nParTypes; ++ip){

      std::string key = std::string("eabs/")+particleName(ip);
      double eabs;
      if(matSection.read(key,eabs) == INTDATA_SUCCESS){
	mat.EABS[ip] = eabs;
	//Update minimum energy eabs for this particle
	if(eabs < minSimEabs[ip])
	  minSimEabs[ip] = eabs;
      }else{
	mat.EABS[ip] = defEabs;	
      }
    }

    // ** C1 and C2
    double C1,C2;
    if(matSection.read("C1",C1) == INTDATA_SUCCESS){
      mat.C1 = C1;
    }
    if(matSection.read("C2",C2) == INTDATA_SUCCESS){
      mat.C2 = C2;
    }

    // ** WCC and WCR
    
    //Set default values
    mat.WCC = std::min(5e3,mat.EABS[PEN_ELECTRON]/100.0);
    mat.WCR = std::min(5e3,mat.EABS[PEN_PHOTON]/100.0);

    //Read values provided by the user
    double WCC, WCR;
    if(matSection.read("WCC",WCC) == INTDATA_SUCCESS){
      mat.WCC = WCC;
    }
    if(matSection.read("WCR",WCR) == INTDATA_SUCCESS){
      mat.WCR = WCR;
    }

    // ** Material file path
    if(matSection.read("filename",simMatFilePaths[imat]) != INTDATA_SUCCESS){
      printf("Error: Unable to read field 'filename' for material '%s'. "
	     "String expected.\n", matName.c_str());
      return -3;
    }
    

    if(verbose > 0){
      printf(" Material '%s' parameters:\n", matName.c_str());
      printf("      C1 =%11.4E       C2 =%11.4E\n",mat.C1,mat.C2);
      printf("     WCC =%11.4E eV,   WCR =%11.4E eV\n",mat.WCC,(mat.WCR > 10.0E0 ? mat.WCR : 10.0E0));
      printf("  electron EABS: %11.4E eV\n",mat.EABS[PEN_ELECTRON]);
      printf("     gamma EABS: %11.4E eV\n",mat.EABS[PEN_PHOTON]);
      printf("  positron EABS: %11.4E eV\n",mat.EABS[PEN_POSITRON]);

      printf("\n Material filename: '%s'.\n",simMatFilePaths[imat].c_str());
      printf("\n Material number: %d.\n",index);
    }
    
  }

  //Configure dummy detector material
  const unsigned detMat = filterMaterialsNames.size()+1;
  pen_material& dummyMat = contextSim.getBaseMaterial(detMat-1);
    
  //Set default C1, C2, WCC and WCR
  dummyMat.C1=0.05;
  dummyMat.C2=0.05;
  dummyMat.WCC=1.0e35;
  dummyMat.WCR=1.0e35;

  //Read absorption energies for each particle type
  for(unsigned ip = 0; ip < constants::nParTypes; ++ip){
      dummyMat.EABS[ip] = 1.0e35;	
  }

  //Init dummy mat with the first material data
  simMatFilePaths[detMat-1].assign(simMatFilePaths[0]);

  
  // ** Filters 
  //**************
  std::vector<std::string> filtersNames;
  std::vector<unsigned> filtersMat;
  err = config.ls("filters",filtersNames);
  if(err != INTDATA_SUCCESS){
    printf("Error: Section 'filters' not found in configuration file.\n");
    return -3;
  }

  if(filtersNames.size() < 1){
    printf("Error: No filters defined in 'filters' section.\n");
    return -3;
  }

  //Resize filter vector information
  filtersMat.resize(filtersNames.size());

  std::vector<bool> usedPos(filtersNames.size());
  std::fill(usedPos.begin(), usedPos.end(), false);
  
  for(const std::string& filterName : filtersNames){

    // Get filter path
    std::string filterPath = std::string("filters/") + filterName;
    pen_parserSection filterSection;
    err = config.readSubsection(filterPath, filterSection);
    if(err != INTDATA_SUCCESS){
      printf("Error: Unable to read filter section '%s'.\n", filterPath.c_str());
      return -3;
    }

    // ** Mat index
    int index;
    err = filterSection.read("material",index);
    if(err != INTDATA_SUCCESS){
      printf("Error: Unable to read 'material' field in filter '%s'\n",
	     filterName.c_str());
      return -3;
    }

    if(index < 0){
      printf("Error on material index assigned to filter '%s' (%d). "
	     "Must be positive\n",filterName.c_str(), index);
    }

    // ** Position
    int position;
    err = filterSection.read("position",position);
    if(err != INTDATA_SUCCESS){
      printf("Error: Unable to read 'position' field in filter '%s'. "
	     "Integer expected\n",
	     filterName.c_str());
      return -3;
    }

    if(position < 0){
      printf("Error on position assigned to filter '%s' (%d). "
	     "Must be positive or zero\n",filterName.c_str(), position);
      return -3;
    }

    if(usedPos[position]){
      printf("Error on position assigned to filter '%s' (%d). "
	     "It is already used\n",filterName.c_str(), position);
      return -3;
    }

    filtersMat[position] = static_cast<unsigned>(index);
    usedPos[position] = true;

    if(verbose > 0){
      printf("\n Filter '%s' parameters:\n", filterName.c_str());
      printf("  Material number: %d.\n",index);
      printf("  Filter position: %d.\n\n",position);
    }
    
  }

  // ** Configure simulation context
  //*********************************
  FILE* fcontext = fopen("context-sim.rep", "w");
  if(contextSim.init(Emax,fcontext,verbose,simMatFilePaths) != PEN_SUCCESS){
    fclose(fcontext);
    printf("Error at simulation context initialization. See context report.\n");
    return -4;
  }
  fclose(fcontext);
  

  // ** Create objective context  
  //*******************************

  //Create a context
  pen_context contextObj;
  
  
  // ** Objective mu values
  //************************
  std::vector<double> objectiveMuRho;
  std::string objectiveMatFilePaths[constants::MAXMAT];  
  std::vector<std::string> objectiveMaterials;
  err = config.ls("objective",objectiveMaterials);
  if(err != INTDATA_SUCCESS){
    printf("Error: Section 'objective' not found in configuration file.\n");
    return -5;
  }

  if(objectiveMaterials.size() < 1){
    printf("Error: No objective values defined in 'objective' section.\n");
    return -5;
  }

  if(objectiveMaterials.size() > constants::MAXMAT){
    printf("Error: Maximum number of objective materials exceeded.\n"
	   "   Maximum: %lu\n"
	   "   Used   : %lu\n",
	   static_cast<unsigned long>(constants::MAXMAT),
	   static_cast<unsigned long>(objectiveMaterials.size()));
    return -5;
  }

  //Create materials in objective context
  const unsigned nObjMats = objectiveMaterials.size();
  errmat = contextObj.setMats<pen_material>(nObjMats);
  if(errmat != 0){
    printf("Error: Unable to create simulation context materials: %d.\n",errmat);
    return -5;
  }
  
  objectiveMuRho.resize(nObjMats);
  unsigned nextMu = 0;
  for(const std::string& matName : objectiveMaterials){

    //Get material path
    std::string objectivePath = std::string("objective/") + matName;
    pen_parserSection objectiveSection;
    err = config.readSubsection(objectivePath, objectiveSection);
    if(err != INTDATA_SUCCESS){
      printf("Error: Unable to read objective material section '%s'.\n", objectivePath.c_str());
      return -5;
    }
    
    // Get mu/rho
    double muRho;
    if(objectiveSection.read("mu",muRho) != INTDATA_SUCCESS){
      printf("Error: Unable to read mu/rho value (field 'mu') from objective "
	     "material '%s'. Double expected.\n",
	     matName.c_str());
      return -5;
    }
    if(muRho <= 0.0){
      printf("Error: Invalid mu/rho value (%E) assigned to objective "
	     "material '%s' in the configuration file.\n",
	     muRho, matName.c_str());
      return -5;
    }
    objectiveMuRho[nextMu] = muRho;

    if(objectiveSection.read("filename",objectiveMatFilePaths[nextMu]) != INTDATA_SUCCESS){
      printf("Error: Unable to read 'filename' from objective "
	     "material '%s'. String expected.\n",
	     matName.c_str());
      return -5;      
    }

    //Configure material parameters
    pen_material& mat = contextObj.getBaseMaterial(nextMu++);
    
    //Set default C1, C2, WCC and WCR
    mat.C1=0.05;
    mat.C2=0.05;
    mat.WCC=defEabs;
    mat.WCR=defEabs;

    for(unsigned ip = 0; ip < constants::nParTypes; ++ip){
      mat.EABS[ip] = minSimEabs[ip];
    }
  }
  
  // ** Configure objective context
  //*********************************
  fcontext = fopen("context-obj.rep", "w");
  if(contextObj.init(Emax,fcontext,verbose,objectiveMatFilePaths) != PEN_SUCCESS){
    fclose(fcontext);
    printf("Error at objective context initialization. See context report.\n");
    return -6;
  }
  fclose(fcontext);

  // ** Number of inspectors
  //****************************
  std::vector<spectrumInspector> inspectors;
  std::vector<FILE*> inspectorsFiles;
  int nInspectors;
  err = config.read("simulation/nInspectors", nInspectors);
  if(err != INTDATA_SUCCESS){
    nInspectors = 10;
  }else{
    if(nInspectors <= 0)
      nInspectors = 10;
  }
  
  inspectors.resize(nInspectors);
  inspectorsFiles.resize(nInspectors);

  //Init inspectors spectrum and parameters
  for(spectrumInspector& inspector : inspectors){
    inspector.spectrum.resize(outputBins, outputEmin, outputDE);
    //Add objective materials
    for(size_t i = 0; i < objectiveMaterials.size(); ++i){
      inspector.addObjectiveMat(i,
				objectiveMaterials[i],
				objectiveMuRho[i],
				contextObj.readBaseMaterial(i).meanZ());
    }
  }

  //-----------------Configuraiton end-------------------------//

  spectrumInspector globalBestResult;
  double globalBestMaxDiff = 1.0e35;

  for(unsigned nFilters = minFilters; nFilters <= filtersNames.size(); ++nFilters){
  
    spectrumInspector bestResult;
    double bestDiffs2 = 1.0e35;

    // ** Init inspectors
    for(size_t iInspector = 0; iInspector < inspectorsFiles.size(); ++iInspector){
      inspectorsFiles[iInspector] = nullptr;
      std::string filename = std::string("filters_") + std::to_string(nFilters) + std::string("_inspector_") +
	std::to_string(iInspector) +  std::string(".dat");
      inspectorsFiles[iInspector] = fopen(filename.c_str(),"w");

      if(inspectorsFiles[iInspector] == nullptr){
	printf("Error: Unable to create file '%s'\n", filename.c_str());
	return -7;
      }

      //Clear previous data
      inspectors[iInspector].clearFilters();      
      
      //Assign first filter sizes randomly
      for(size_t iFilter = 0; iFilter < nFilters; ++iFilter){
	const unsigned index = (nFilters-1)-iFilter;
	const double thickness = (0.05 + rand01(gen)* 0.25)/static_cast<double>(nFilters);

	const unsigned imat = filtersMat[index];

	if(imat == 0){
	  inspectors[iInspector].
	    filters.emplace_back(thickness,
				 0.0,
				 0);
	}else{
	  inspectors[iInspector].
	    filters.emplace_back(thickness,
				 contextSim.getBaseMaterial(imat-1).meanZ(),
				 imat);
	}
      }

      
      inspectors[iInspector].printRowHeader(inspectorsFiles[iInspector]);
    }
  
    for(size_t iIter = 0; iIter < nIter; ++iIter){

      if(verbose > 1){
	printf("\n **** Iteration %d\n\n", static_cast<unsigned>(iIter));
      }
    
      for(int iInspector = 0; iInspector <  static_cast<int>(inspectors.size()); ++iInspector){

	spectrumInspector& inspector = inspectors[iInspector];
      
	// ** Configure geometry
	//***********************
	pen_parserSection geoConfig;

	//Find required energy to cross the whole filters
	//except the first one
	std::vector<double> filterPhotonRangeEnergy(inspector.filters.size());
	filterPhotonRangeEnergy[0] = 0.0;
	for(size_t ifilter = 1; ifilter < inspector.filters.size(); ++ifilter){

	  //Get filter pair information
	  const filter& f = inspector.filters[ifilter];

	  //Avoid void filters
	  if(f.imat == 0){
	    filterPhotonRangeEnergy[ifilter] = 0.0;
	    continue;
	  }

	  //Calculate the required gamma range to absorb 99% of incident photons
	  //
	  // I = I0 exp(-l * mu)
	  // ln(I/I0) = -l * mu
	  // mu = - ln(I/I0) / l
	  // range = 1/mu
	  //
	  const double objLogTrans = -log(0.01);
	  const double objMu = objLogTrans / f.thickness;
	  const double objRange = 1.0/objMu;
	  
	  double lowE = contextSim.getMatEABS(f.imat-1, PEN_PHOTON);
	  double topE = Emax;
	  contextSim.findRange(lowE, topE, objRange,
			       PEN_PHOTON, f.imat-1);
	  if(topE == Emax)
	    filterPhotonRangeEnergy[ifilter] = 1.0e35;
	  else
	    filterPhotonRangeEnergy[ifilter] = lowE;
	  
	}
	
	//Construct the geometry according to inspector filters specifications
	const unsigned nSubFilters = 10;
	for(size_t ifilter = 0; ifilter < inspector.filters.size(); ++ifilter){

	  //Get filter pair information
	  const filter& f = inspector.filters[ifilter];
	  
	  //Iterate over subfilters
	  const double subFilterSize =
	    f.thickness / static_cast<double>(nSubFilters);
	  for(unsigned iSubFilter = 0; iSubFilter < nSubFilters; ++iSubFilter){

	    const std::string basePath =
	      std::string("filters/") + std::to_string(ifilter) +
	      std::string("_") + std::to_string(iSubFilter);

	    
	    std::string auxStr = basePath + std::string("/position");
	    int pos = static_cast<int>(ifilter*nSubFilters + iSubFilter);
	    geoConfig.set(auxStr, pos);
	    auxStr = basePath + std::string("/width");
	    geoConfig.set(auxStr, subFilterSize);
	    auxStr = basePath + std::string("/material");
	    geoConfig.set(auxStr, static_cast<int>(f.imat));


	    if(f.imat > 0){
	      double eabs = defEabs;
	      //If this subfilter has subfilters beyond, recalculate eabs
	      if(iSubFilter < nSubFilters -1){
		//Calculate the required gamma range to absorb 99.99% of created photons
		//in a length equal to the remaining subfilters thickness 
		//
		const double objLogTrans = -log(0.0001);
		const double length =
		  subFilterSize*static_cast<double>(nSubFilters - (iSubFilter+1));
		const double objMu = objLogTrans / length;
		const double objRange = 1.0/objMu;

		//Calculate and set absorption energy for photons according to filter thickness
		double lowE = contextSim.getMatEABS(f.imat-1, PEN_PHOTON);
		double topE = Emax;
	  
		contextSim.findRange(lowE, topE, objRange,
				     PEN_PHOTON, f.imat-1);
		if(topE == Emax)
		  eabs = 1.0e35;
		else
		  eabs = lowE;
	      }

	      //Check if generated photons can cross next filters
	      if(ifilter < inspector.filters.size()-1){
		double nextEabs =
		  *std::max_element(filterPhotonRangeEnergy.begin() + ifilter + 1,
				    filterPhotonRangeEnergy.end());
		if(eabs < nextEabs)
		  eabs = nextEabs;
	      }

	      //Set absorption energies in configuration
	      std::string baseEabsPath = basePath + std::string("/eabs/");
	      for(unsigned ip = 0; ip < constants::nParTypes; ++ip){

		//Set local absorption energy
		auxStr = baseEabsPath + particleName(ip);
		geoConfig.set(auxStr, eabs);
	  
	      }
	    }
	    
	  }
	
	}

	//Add a final filter as detector
	geoConfig.set("filters/detector/position",
		      static_cast<int>(nSubFilters*inspector.filters.size()));
	geoConfig.set("filters/detector/width", 1.0);
	geoConfig.set("filters/detector/material", static_cast<int>(detMat));      
      
	//Configure the geometry
	pen_filterGeo geometry;
	err = geometry.configure(geoConfig,1);
	if(err != 0){
	  printf("Error configuring geometry.\n");
	  return -3;
	}

	//Set the geometry to the simulation context
	contextSim.setGeometry(&geometry);

	//Update eabs information
	contextSim.updateEABS();

	if(verbose > 2){
	  printf("\n Geometry config: \n%s\n\n", geoConfig.stringify().c_str());
	}
            
	// ** Start the simulation  
	//**************************

	//Adjust the tolerance according to the difference between
	//objective measurment and the obtained in last iteration
	double objectiveTol = 100.0*std::fabs(inspector.relMaxDiff())/2.0;
	if(objectiveTol > 100.0)
	  objectiveTol = 30.0;
	else if(objectiveTol < tolerance)
	  objectiveTol = tolerance;

	if(verbose > 1)
	  printf(" Objective tolerance for inspector %d spectrum: %.3f%%\n",
		 iInspector,objectiveTol);
	
#ifdef _PEN_USE_THREADS_
	std::vector<std::thread> threads;	
	const double localTol = objectiveTol*sqrt(nThreads);
	for(size_t ith = 0; ith < nThreads; ++ith){
	
	  unsigned long long localHists = nHists/nThreads;
	  if(ith == 0){
	    localHists += nHists % nThreads;
	  }
	
	  threads.push_back(std::thread(simulate,
					localHists,
					localTol,
					&inspector.system,
					&contextSim,
					&psfs[ith],
					detMat,
					&results[ith],
					&seeds1[ith],&seeds2[ith]));
	
	}

	//Wait threads finish
	for(std::thread& t : threads){
	  t.join();
	}

	//Clear thread vector
	threads.clear();

	for(size_t ith = 1; ith < nThreads; ++ith){
	  results[0].sum(results[ith]);
	}
#else
	simulate(nHists,objectiveTol,
		 &inspector.system,
		 &contextSim,&psfs[0],
		 detMat,&globResults,
		 seeds1[0],seeds2[0]);      
#endif

	//Save resulting spectrum with errors
      
	std::pair<std::vector<double>, std::vector<double>> res = globResults.results();
      
	inspector.spectrum.weight.swap(res.first);
	inspector.spectrum.eWeight.swap(res.second);
      
	//Calculate cummulative function of resulting spectrum
	if(inspector.spectrum.updateCummulative() != 0){
	  
	  printf("Empty spectrum. Resize filters and repeat the simulation.\n\n");

	  for(filter& f : inspector.filters){
	    f.thickness /= 2.0;
	  }
	  --iInspector; 	  
	  continue;
	}
      
	// ** Mu Calculus
	//*****************
      
	//Calculate mu values for the resulting spectrum

#ifdef _PEN_USE_THREADS_

	for(size_t ith = 0; ith < nThreads; ++ith){
	
	  unsigned long long localHists = nHists/nThreads;
	  if(ith == 0){
	    localHists += nHists % nThreads;
	  }
	
	  threads.push_back(std::thread(calculateMu,
					localHists,
					localTol,
					nObjMats,
					&inspector.spectrum,
					&contextObj,
					&muThreads[ith],
					&seeds1[ith], &seeds2[ith]));
	}

	//Wait until all threads finish
	for(std::thread& t : threads){
	  t.join();
	}

	//Clear thread vector
	threads.clear();      

	//Sum results
	for(size_t ith = 1; ith < nThreads; ++ith){
	  muThreads[0].sum(muThreads[ith], nObjMats);
	}
      
#else
	calculateMu(nHists,
		    tolerance,
		    nObjMats,
		    &inspector.spectrum,
		    &contextObj,
		    &globMuResults,
		    &seeds1[0], &seeds2[0]);
#endif

	//Calculate final mu and differences values for each material
	const double iNIter = 1.0/static_cast<double>(globMuResults.nIter);
	for(unsigned im = 0; im < nObjMats; ++im){

	  const double meanRange = globMuResults.RANGET[im]*iNIter;
	  const double range2 = meanRange*meanRange;
	  const double meanRange2 = globMuResults.RANGET2[im]*iNIter;
	  double eRange = (meanRange2 - range2)*iNIter;
	  if(eRange > 0)
	    eRange = sqrt(eRange);
	  else
	    eRange = 0.0;

	  const pen_material& mat = contextObj.readBaseMaterial(im);
	  const double q = 1.0/(mat.RHO*meanRange);
	  double eq = eRange/(mat.RHO*range2);
	    
	  inspector.setMu(im,q,eq);
	}

	inspector.printRow(inspectorsFiles[iInspector]);

	//Check the finish condition
	if(100.0*fabs(inspector.relMaxDiff()) <= tolerance){
	  if(verbose > 1){
	    printf("\nReached objective tolerance: %.3f%%\n",
		   100.0*fabs(inspector.relMaxDiff()));
	    inspector.print(true,true,false);	    
	  }
	  
	  if(objectiveTol > tolerance){
	    //Repeat the simulation with better tolerance with no geometry changes
	    if(verbose > 1){
	      printf("\nThe spectrum has been obtained with low precission (%.3f%%).\n"
		     " The simulations will be repeated to obtain the specified precission: %.3f%%\n",
		     objectiveTol,tolerance);
	    }
	    --iInspector;
	    continue;
	  }
	  break;
	}
      }

      //Find the inspector with greater agreement
      unsigned best = 0;
      double minSumDiff2 = inspectors[0].sumDiffs2();
      for(size_t i = 1; i < inspectors.size(); ++i){
	const double sumDiff2 = inspectors[i].sumDiffs2();
    
	if(sumDiff2 < minSumDiff2){
	  minSumDiff2 = sumDiff2;
	  best = i;
	}
      
      }

      //Check if the best inspector in this iteration is better than the old one
      if(minSumDiff2 < bestDiffs2){
	bestDiffs2 = minSumDiff2;
	bestResult = inspectors[best];
      }

      if(verbose > 1){
	printf("\n\n ** Best iteration result:\n\n");
	bestResult.print(true,true,false);
      }      

      //Check the finish condition
      if(100.0*fabs(bestResult.relMaxDiff()) <= tolerance){
	break;
      }

      //Select the worst mu value to correct it in each inspector
      pen_rand random;
      random.setSeeds(seeds1[0], seeds2[0]);
      for(spectrumInspector& inspector : inspectors){
	double maxDiff = inspector.relMaxDiff();

	//Select the filter to modify
	unsigned ifilter =
	  static_cast<unsigned>(random.rand()*static_cast<double>(inspector.nFilters()));
	ifilter %= inspector.nFilters(); //Just in case

	//If is a void material, change the resize direction
	if(inspector.filters[ifilter].imat == 0)
	  maxDiff = -maxDiff;

	//Calculate change factor and limit it
	double changeFactor = maxDiff*( 2.0*random.rand() );
	if(fabs(changeFactor) > 0.5){
	  if(changeFactor < 0.0)
	    changeFactor = -0.5;
	  else
	    changeFactor = 0.5;
	}

	inspector.filters[ifilter].thickness *= 1.0 + changeFactor;
	if(inspector.filters[ifilter].thickness < 1.0e-5)
	  inspector.filters[ifilter].thickness = 1.0e-5;
	
      }
      random.getSeeds(seeds1[0], seeds2[0]);
    }

    for(FILE* inspectorFile : inspectorsFiles){
      fclose(inspectorFile);
    }

    std::string resFilename = std::string("best_") +
      std::to_string(nFilters) +
      std::string("_filters.dat");
    
    FILE* fbest = fopen(resFilename.c_str(), "w");
    bestResult.fprint(fbest);
    fclose(fbest);

    if(bestResult.relMaxDiff() < globalBestMaxDiff){
      globalBestMaxDiff = bestResult.relMaxDiff();
      globalBestResult = bestResult;
    }
  }
  
  FILE* fbest = fopen("bestResult.dat", "w");
  globalBestResult.fprint(fbest);
  fclose(fbest);

  return 0;
}

void simulate(const unsigned long long nIter,
	      const double tolerance,
	      const systemGeo* system,
	      const pen_context* pSimContext,
	      psfMemory_specificSampler* sampler,
	      const unsigned detMat,
	      measurmentRegular<double>* results,
	      int* seed1, int* seed2){

  //Reset results
  results->reset();

  //Reset sampler
  sampler->reset();
  
  //Create random generator
  pen_rand random;
  random.setSeeds(*seed1,*seed2);
  
  const pen_context& context = *pSimContext;

  //Create particle stacks for electrons, positrons and gamma
  pen_particleStack<pen_particleState> stackE;
  pen_particleStack<pen_particleState> stackP;
  pen_particleStack<pen_state_gPol> stackG;
  

  //Create particle instances
  pen_betaE betaE(context,stackE,stackG);
  pen_gamma gamma(context,stackE,stackP,stackG);
  pen_betaP betaP(context,stackE,stackG,stackP);

  pen_particleState genState;
  unsigned long long hist = 0;
  
  while(hist < nIter){

    pen_KPAR genKpar;
    unsigned long long dhist;
    sampler->sample(genState,genKpar,dhist,random);
    if(genKpar == pen_KPAR::ALWAYS_AT_END){
      //End of psf
      break;
    }

    //Add source offset
    genState.Z += system->distPSF2Filter;

    //Increase history counter
    hist += dhist;

    if(dhist > 0){
      //End counting for this history
      results->checkPoint();      
    }

    //Check if it is a gamma
    if(genKpar != pen_KPAR::PEN_PHOTON){
      //Skip it
      continue;
    }
    
    context.readGeometry()->locate(genState);
    
    //Copy generated state
    stateCopy(gamma.getState(),genState);

    //Init Page variables
    gamma.page0();
    
    //Try to move generated particle to geometry system
    if(move2geo(gamma)){

      //Simulate sampled particle
      simulatePart(gamma,random);

      //Check if the particle has escaped the geometry in -Z direction
      const pen_particleState& state = gamma.readState();
      
      if(state.MAT == detMat){
	//Count this particle in the final spectrum
	//results->add(state.E, state.WGHT);

	//Check if the particle can be detected
	if(system->isDetected(state)){
	  //Count this particle in the final spectrum
	  results->add(state.E, state.WGHT);
	}	
      }
    }
    
    //Secondary particles
    for(;;){

      //Get the state and number of stacked particles
      unsigned nBetaE = stackE.getNSec();
      unsigned nGamma = stackG.getNSec();
      unsigned nBetaP = stackP.getNSec();

      unsigned nBetaE05 = 0;
      unsigned nGamma05 = 0;
      unsigned nBetaP05 = 0;
      
      //Check end of simulation
      bool remaining = false;
      if(nBetaE > 0){
	remaining = true;
	nBetaE05 = nBetaE/2+1;
      }
      if(nGamma > 0){
	remaining = true;
	nGamma05 = nGamma/2+1;
      }
      if(nBetaP > 0){
	remaining = true;
	nBetaP05 = nBetaP/2+1;
      }

      //Check if the stacks are empty
      if(!remaining)
	break;
      
      
      //Simulate half stacked particles
      unsigned nbetaEsim = 0;
      while(nbetaEsim < nBetaE05){
	stackE.get(betaE.getState());

	//Update body and material
	betaE.updateBody();
	betaE.updateMat();
	  	  
	//Simulate particle
	simulatePart(betaE,random);	
	nbetaEsim++;
      }
	
      unsigned ngammasim = 0;
      while(ngammasim < nGamma05){
	stackG.get(gamma.getState());

	//Update body and material
	gamma.updateBody();
	gamma.updateMat();
	    
	//Simulate particle
	simulatePart(gamma,random);
	ngammasim++;

	const pen_particleState& state = gamma.getState();
	if(state.MAT == detMat){

	  //Count this particle in the final spectrum
	  //results->add(state.E, state.WGHT);
	  
	  //Check if the particle can be detected
	  if(system->isDetected(state)){	  
	    //Count this particle in the final spectrum
	    results->add(state.E, state.WGHT);
	  }

	}
      }
	
      unsigned nbetaPsim = 0;
      while(nbetaPsim < nBetaP05){
	stackP.get(betaP.getState());

	//Update body and material
	betaP.updateBody();
	betaP.updateMat();
	  
	//Simulate particle
	simulatePart(betaP,random);	
	nbetaPsim++;
      }
    }

    //Check if the simulation can finish
    if(hist > 10000 && hist % 1000 == 0){
      double maxAbsErel = results->maxAbsErel(0.95);
      if(maxAbsErel <= tolerance/100.0)
	break;
    }
  }

  //Update seeds
  random.getSeeds(*seed1,*seed2);
}


void calculateMu(const unsigned long long nIter,
		 const double tolerance,
		 const unsigned nMat,
		 const spectrumStruct* spectrum,
		 pen_context* contextObj,
		 calcMuStruct* results,
		 int* seed1, int* seed2){

  //Init random generator
  pen_rand random;
  random.setSeeds(*seed1,*seed2);

  //Reset results structure
  results->reset(nMat);

  results->nIter = nIter;
  
  for(unsigned long long i = 1; i <= nIter; ++i){
	
    const double E0 = spectrum->sample(random);
	
    //Iterate over materials
    for(unsigned im = 0; im < nMat; ++im){

      //Calculate and save range contribution
      const double RANGE = contextObj->range(E0,PEN_PHOTON,im);

      results->RANGET[im]  += RANGE;
      results->RANGET2[im] += RANGE*RANGE;
    }

    //Check error
    if(i % 10000 == 0 || i == nIter){
      double maxErel = 0.0;
      const double iNIter = 1.0/static_cast<double>(i);
      for(size_t im = 0; im < nMat; ++im){

	const double meanRange = results->RANGET[im]*iNIter;
	const double range2 = meanRange*meanRange;
	const double meanRange2 = results->RANGET2[im]*iNIter;
	double eRange = (meanRange2 - range2)*iNIter;
	if(eRange > 0)
	  eRange = sqrt(eRange);
	else
	  eRange = 0.0;

	const pen_material& mat = contextObj->readBaseMaterial(im);
	const double q = 1.0/(mat.RHO*meanRange);
	double eq = eRange/(mat.RHO*range2);
	    
	//inspector->muRhos[im] = q;
	//inspector->eMuRhos[im] = 100.0*eq/q;

	double eRel = 100.0*eq/q;
	if(maxErel < eRel)
	  maxErel = eRel;
      }

	
      if(maxErel <= tolerance){
	results->nIter = i;
	break;
      }
    }
  }

  //Update seeds
  random.getSeeds(*seed1,*seed2);
  
}
