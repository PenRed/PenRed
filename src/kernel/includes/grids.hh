#ifndef _PEN_MATH_GRIDS_
#define _PEN_MATH_GRIDS_

#include <array>
#include <vector>
#include <utility>
#include <algorithm>
#include <numeric>
#include <functional>
#include <cstdio>
#include <istream>
#include <string>
#include <limits>
#include <sstream>
#include <cmath>

namespace penred {

  inline std::string triml(const std::string& strin){
    const std::string delimiters = " \n\r\t\f\v";
    size_t first = strin.find_first_not_of(delimiters);
    return (first == std::string::npos) ? "" : strin.substr(first);
  }

  inline std::string trimr(const std::string& strin){
    const std::string delimiters = " \n\r\t\f\v";
    size_t last = strin.find_last_not_of(delimiters);
    return (last == std::string::npos) ? "" : strin.substr(0,last+1);
  }

  inline std::string trim(const std::string& strin){
    return trimr(triml(strin));
  }
  
  namespace grids {

    enum errors{
      SUCCESS = 0,
      INVALID_LIMITS,
      LIMITS_NOT_FOUND,
      INVALID_NUMBER_OF_BINS,
      DIMENSION_MISMATCH,
      DIMENSION_OUT_OF_RANGE,
      DIMENSION_REPEATED,
      DIMENSION_NOT_FOUND,
      UNEXPECTED_DIMENSIONS,
      MIN_BIN_GREATER_THAN_MAX,
      BIN_OUT_OF_RANGE,
      NUMBER_OF_BINS_NOT_FOUND,
    };

    constexpr const char* errorToString( const unsigned i ){
      switch(i){
      case SUCCESS: return "Success";
      case DIMENSION_MISMATCH: return "Dimensions mismatch";
      case DIMENSION_OUT_OF_RANGE: return "Dimension out of range";
      case DIMENSION_REPEATED: return "Dimension repeated";
      case DIMENSION_NOT_FOUND: return "Dimensions not found";
      case UNEXPECTED_DIMENSIONS: return "Unexpected dimensions";
      case INVALID_LIMITS: return "Invalid limits";
      case LIMITS_NOT_FOUND: return "Limits not found";
      case INVALID_NUMBER_OF_BINS: return "Invalid number of bins";
      case MIN_BIN_GREATER_THAN_MAX: return "Minimum bin grater than maximum";
      case BIN_OUT_OF_RANGE: return "Bin out of range";
      case NUMBER_OF_BINS_NOT_FOUND: return "Number of bins not found";
      default: return "Unknown error";
      };
    }
    
    // ============================================================================
    //  Grid contract (enforced by convention)
    // ============================================================================
    //
    //  A grid type G used as the `Grid` template parameter of
    //  multiDimension<dim,G>, results<T,dim,G>, measurement<T,dim,G> must provide:
    //
    //  --- Index geometry (normally inherited from BinningBase<dim>) -------------
    //      std::array<unsigned long, dim> nBins;
    //      std::array<unsigned long, dim> binsPerIncrement;
    //      unsigned long                  totalBins;
    //      unsigned long globalIndex(const std::array<unsigned long,dim>&) const;
    //      bool sameBins(const BinningBase& o) const;
    //      unsigned long effectiveDim() const;
    //      template<size_t profDim> getProfileBins(const std::array<size_t, profDim>& profDimsIndex,
    //                                  const std::array<std::pair<unsigned long, unsigned long>, dim>& binLimits,
    //                                  std::array<unsigned long, profDim>& profileBins) const
    //
    //  --- Coordinate geometry (each grid implements these) ----------------------
    //      std::array<std::pair<double,double>, dim> limits;
    //      std::array<unsigned long, dim>
    //          findBin(const std::array<double,dim>& pos) const;
    //          // Bin index containing pos. Must clamp to [0, nBins[d]-1].
    //          // Out-of-range rejection is the CALLER's job.
    //
    //
    //      double binLower(size_t d, unsigned long i) const;  // left edge of bin i
    //      double binUpper(size_t d, unsigned long i) const;  // right edge of bin i
    //      double lowEdge (size_t d) const;                   // = binLower(d,0)
    //      double highEdge(size_t d) const;                   // = binUpper(d,nBins-1)
    //      double binWidth(size_t d, unsigned long i) const   // Returns bin's width
    //
    //  --- bracket(pos, lo, hi, frac)
    //  For each dimension d:
    //    lo[d]   = bin index containing pos[d], clamped to [0, nBins[d]-1]
    //    hi[d]   = lo[d] + 1 if lo[d] < nBins[d]-1, else lo[d]  (no right neighbour)
    //    frac[d] = position of pos[d] inside [binLower(lo[d]), binUpper(lo[d])],
    //              in [0,1].  frac = 0 at the left edge, 1 at the right edge.
    //
    //  When hi[d] == lo[d] (last bin or out of range), the interpolated value is constant
    //  data[lo[d]] regardless of frac[d]; frac[d] is still reported truthfully
    //  (1.0 exactly at highEdge) for callers that inspect it directly.
    //      void bracket(const std::array<double,dim>& pos,
    //                   std::array<unsigned long,dim>& lower,
    //                   std::array<unsigned long,dim>& upper,
    //                   std::array<double,dim>& frac) const;
    //          // lower[d] = bin containing pos[d];
    //          // upper[d] = min(lower[d]+1, nBins[d]-1);
    //          // frac[d]  in [0,1] = position of pos[d] inside [lower,upper];
    //          // if lower==upper, frac must be 0.
    //
    //  --- Generate initialization edges
    //      std::array<std::vector<double>, dim> getEdges(const bool last = true) const
    //      std::array<std::vector<double>, dim> getInitEdges() const
    //            template<size_t profDim> int getProfileInit(const std::array<size_t, profDim>& profDimsIndex,
    //                     const std::array<std::pair<unsigned long, unsigned long>, dim>& binLimits,
    //                     std::array<unsigned long, profDim>& profileBins,
    //                     std::array<std::vector<double>, profDim>& profEdges) const
    //
    //  For RegularGrid: profEdges[i] has size 2 (the {low, high} limits
    //  of the sub-range), matching RegularGrid::init's expectation.
    //  For VariableGrid: profEdges[i] has size profileBins[i]+1 (the full
    //  edge list), matching VariableGrid::init's expectation.    
    //
    //  --- Identity (used by subtract/relative compatibility checks) -------------
    //      bool sameAs(const G& other) const;
    //          // true iff both grids describe the SAME binning
    //          // (same nBins AND same edges/limits).
    //
    //  --- Binning I/O -----------------------------------------------------------
    //      void printDims(FILE* fout, bool printOnlyEffective) const;
    //      int  loadDims(std::istream& in);
    //          // Must set nBins and the grid's own representation, then call
    //          // this->computeIndexGeometry() before returning errors::SUCCESS.
    //
    //  --- Dimension rebinding (used by profileByBins / to1D / to2D) -----------
    //      template<size_t newDim> using rebind = <SameGridTemplate><newDim>;
    //          // Must map the grid's family to the same family at newDim:
    //          //   RegularGrid<dim>::rebind<N>  == RegularGrid<N>
    //          //   VariableGrid<dim>::rebind<N> == VariableGrid<N>
    //          // The measurements layer uses grids::rebind_grid_t<profDim, Grid>
    //          // to declare the profile output type.
    //
    //  Reference implementations: RegularGrid<dim>, VariableGrid<dim>.
    // ============================================================================

    template<size_t newDim, class Grid>
    using rebind_t = typename Grid::template rebind<newDim>;
    
    template<size_t dim>
    struct BinningBase {
      std::array<unsigned long, dim> nBins{};
      std::array<unsigned long, dim> binsPerIncrement{};
      unsigned long totalBins = 0;

      void computeIndexGeometry() {
        totalBins = 1ul;
        for(size_t i = 0; i < dim; ++i) totalBins *= nBins[i];
        binsPerIncrement[0] = 1;
        for(size_t i = 1; i < dim; ++i)
          binsPerIncrement[i] = binsPerIncrement[i-1] * nBins[i-1];
      }

      void clearBins(){
        for(unsigned long& b : nBins)
          b = 1;
        computeIndexGeometry();
      }      

      unsigned long globalIndex(const std::array<unsigned long, dim>& idx) const {
        unsigned long g = idx[0];
        for(size_t i = 1; i < dim; ++i) g += idx[i] * binsPerIncrement[i];
        return g;
      }

      bool sameBins(const BinningBase& o) const {
        for(size_t d = 0; d < dim; ++d){
          if(nBins[d] != o.nBins[d]) return false;
        }
        return true;
      }

      inline unsigned long effectiveDim() const {
        //Calculate effective dimensions
        unsigned long effectiveDims = dim;
        for(size_t i = 0; i < dim; ++i){
          if(nBins[i] == 1){
            --effectiveDims;
          }
        }
        return effectiveDims;
      }

      template<size_t profDim>
      int getProfileBins(const std::array<size_t, profDim>& profDimsIndex,
                         const std::array<std::pair<unsigned long, unsigned long>, dim>& binLimits,
                         std::array<unsigned long, profDim>& profileBins) const {
        //Returns the number of bins for a profile with the specified bin ranges
        //
        // profDimsIndex: Original dimension's indexes used to generate the profile. Dimensions not included will be aggregated
        // binLimits:     Original grid range where the profile is extracted
        // profileBins:   Number of bins in each profile dimension
        //

        //Check dimensions
        if(profDim == 0 || profDim > dim){
          return errors::DIMENSION_OUT_OF_RANGE;
        }        
        for(size_t i = 0; i < profDim; ++i){
          if(profDimsIndex[i] >= dim)
            return errors::DIMENSION_OUT_OF_RANGE;
        }
        for(size_t i = 0; i < profDim-1; ++i){
          for(size_t j = i+1; j < profDim; ++j){
            if(profDimsIndex[i] == profDimsIndex[j])
              return errors::DIMENSION_REPEATED;
          }
        }
        for(size_t i = 0; i < dim; ++i){
          if(binLimits[i].first >= binLimits[i].second)
            return errors::MIN_BIN_GREATER_THAN_MAX;
          if(binLimits[i].second > this->nBins[i])
            return errors::BIN_OUT_OF_RANGE;
        }

        //Fill profile bins
        for(size_t i = 0; i < profDim; ++i){
          profileBins[i] =
            binLimits[profDimsIndex[i]].second - binLimits[profDimsIndex[i]].first;
        }

        return errors::SUCCESS;
      }

      void printBinning(FILE* fout, const bool printOnlyEffective) const {
        if(printOnlyEffective){
          fprintf(fout, "# Dimensions:\n"
                  "# %lu\n"
                  "# Number of bins:\n"
                  "#",
                  effectiveDim());
        }
        else{
          fprintf(fout, "# Dimensions:\n"
                  "# %lu\n"
                  "# Number of bins:\n"
                  "#",
                  static_cast<unsigned long>(dim));
        }
	
        for(size_t i = 0; i < dim; ++i){
          if(printOnlyEffective && nBins[i] == 1){
            continue;
          }
          fprintf(fout, " %lu", nBins[i]);
        }
        fprintf(fout, "\n");        
      }

      int loadBinning(std::istream& in, unsigned long& readDim,
                      std::vector<unsigned long>& readNBins,
                      unsigned long& readEffectiveDim,
                      bool& readEffective){

        //Read dimension
        std::string line;
        std::getline(in, line);
        if(sscanf(line.c_str(), " %*c %lu", &readDim) != 1){
          return errors::DIMENSION_NOT_FOUND;
        }
        
        // ** Number of bins
        bool found = false;
        while(std::getline(in, line)){
          if(line.find("# Number of bins:") == 0){
            found = true;
            break;
          }
        }


        if(!found){
          return errors::NUMBER_OF_BINS_NOT_FOUND;
        }

        //Remove comment character
        in.ignore(std::numeric_limits<std::streamsize>::max(), '#');        

        //Read number of bins for each dimension
        readNBins.resize(readDim);
        for(size_t i = 0; i < readDim; ++i){
          in >> readNBins[i];
          if(!in){
            return errors::NUMBER_OF_BINS_NOT_FOUND;
          }
          if(readNBins[i] == 0){
            return errors::INVALID_NUMBER_OF_BINS;
          }
        }

        //Count effective dimensions
        readEffectiveDim = 0;
        for(size_t i = 0; i < readDim; ++i){
          if(readNBins[i] > 1){
            ++readEffectiveDim;
          }
        }

        // * Check dimensions
        readEffective = false;
        if(readDim > dim){
          if(readEffectiveDim <= dim)
            readEffective = true;
          else
            return errors::DIMENSION_MISMATCH;
        }

        //Set number of bins
        clearBins();
        if(readEffective){
          size_t j = 0;
          for(size_t i = 0; i < readDim; ++i){
            if(readNBins[i] > 1){
              nBins[j++] = readNBins[i];
            }
          }
        }
        else{
          for(size_t i = 0; i < readDim; ++i){
            nBins[i] = readNBins[i];
          }
        }

        //Calculate bin geometry
        computeIndexGeometry();

        //Skip line
        in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        return errors::SUCCESS;
      }
    };

    // ---------------------------------------------------------------------------
    //  RegularGrid: uniform bin width per dimension
    // ---------------------------------------------------------------------------
    template<size_t dim>
    class RegularGrid : public BinningBase<dim> {
    public:
      template<size_t newDim>
      using rebind = RegularGrid<newDim>;
      
      std::array<std::pair<double,double>, dim> limits{};
      std::array<double, dim> binWidths{};

      // Build from number of bins + limits. Returns an int error code
      int init(const std::array<unsigned long, dim>& nBinsIn,
               const std::array<std::vector<double>, dim>& limitsIn) {
        for(size_t d = 0; d < dim; ++d){
          if(limitsIn[d].size() != 2) return errors::UNEXPECTED_DIMENSIONS;
          if(limitsIn[d][0] >= limitsIn[d][1]) return errors::INVALID_LIMITS;
          if(nBinsIn[d] == 0)                         return errors::INVALID_NUMBER_OF_BINS;
        }
        this->nBins  = nBinsIn;
        for(size_t d = 0; d < dim; ++d){
          limits[d].first = limitsIn[d][0];
          limits[d].second = limitsIn[d][1];
          binWidths[d] = (limits[d].second - limits[d].first)
            / static_cast<double>(this->nBins[d]);
        }
        this->computeIndexGeometry();
        return errors::SUCCESS;
      }

      std::array<unsigned long, dim>
      findBin(const std::array<double,dim>& pos) const {
        std::array<unsigned long, dim> idx;
        for(size_t d = 0; d < dim; ++d){
          const double f = (pos[d] - limits[d].first) / binWidths[d];
          if(f < 0.0) idx[d] = 0;
          else if(f >= static_cast<double>(this->nBins[d])) idx[d] = this->nBins[d] - 1;
          else idx[d] = static_cast<unsigned long>(f);
        }
        return idx;
      }

      void bracket(const std::array<double,dim>& pos,
                   std::array<unsigned long,dim>& lo,
                   std::array<unsigned long,dim>& hi,
                   std::array<double,dim>& frac) const {
        for(size_t d = 0; d < dim; ++d){
          double f = (pos[d] - this->limits[d].first) / this->binWidths[d];
          // Clamp to the domain so frac stays in [0,1]
          if(f < 0.0) {
            // Below lower limit, use first bin
            lo[d] = 0;
            hi[d] = 0;
            frac[d] = 0.0;
          }
          else if(f >= static_cast<double>(this->nBins[d])){
            // Above upper limit, use last bin
            lo[d] = this->nBins[d]-1;
            hi[d] = this->nBins[d]-1;
            frac[d] = 1.0;
          }
          else{
            //Within range
            lo[d] = static_cast<unsigned long>(f);
            //if(lo[d] >= this->nBins[d]) lo[d] = this->nBins[d] - 1;

            hi[d] = std::min(lo[d] + 1, this->nBins[d] - 1);
            
            frac[d] = f - static_cast<double>(lo[d]);
          }
        }
      }

      double binLower(size_t d, unsigned long i) const {
        return limits[d].first + static_cast<double>(i) * binWidths[d];
      }
      double binUpper(size_t d, unsigned long i) const {
        return limits[d].first + static_cast<double>(i + 1) * binWidths[d];
      }
      double lowEdge (size_t d) const { return limits[d].first;  }
      double highEdge(size_t d) const { return limits[d].second; }
      double binWidth(size_t d, unsigned long /*i*/) const {
        return binWidths[d];
      }

      std::array<std::vector<double>, dim> getEdges(const bool last = true) const {
        std::array<std::vector<double>, dim> edg;
        for(size_t d = 0; d < dim; ++d){
          edg[d].resize(last ? this->nBins[d]+1 : this->nBins[d]);
          for(size_t i = 0; i < this->nBins[d]; ++i){
            edg[d][i] = binLower(d,i);
          }
          if(last) edg[d][this->nBins[d]] = binLower(d,this->nBins[d]);
        }
        return edg;
      }
      std::array<std::vector<double>, dim> getInitEdges() const {
        std::array<std::vector<double>, dim> edg;
        for(size_t d = 0; d < dim; ++d){
          edg[d].resize(2);
          edg[d][0] = limits[d].first;
          edg[d][1] = limits[d].second;
        }
        return edg;
      }
      
      template<size_t profDim>
      int getProfileInit(const std::array<size_t, profDim>& profDimsIndex,
                         const std::array<std::pair<unsigned long, unsigned long>, dim>& binLimits,
                         std::array<unsigned long, profDim>& profileBins,
                         std::array<std::vector<double>, profDim>& profEdges) const {

        int errBin = this->getProfileBins(profDimsIndex, binLimits, profileBins);
        if(errBin != errors::SUCCESS)
          return errBin;

        //Extract edges
        for(size_t i = 0; i < profDim; ++i){
          //Get dimension index in the original grid
          size_t profIndex = profDimsIndex[i];

          //Save edges
          profEdges[i].resize(2);
          profEdges[i][0] = this->limits[profIndex].first + 
            static_cast<double>(binLimits[profIndex].first)*this->binWidths[profIndex];
          profEdges[i][1] = this->limits[profIndex].first + 
            static_cast<double>(binLimits[profIndex].second)*this->binWidths[profIndex];
        }
        
        return errors::SUCCESS;
      }
      
      const std::array<std::pair<double, double>, dim>& readLimits() const { return limits; }

      bool sameAs(const RegularGrid& o, const double threshold = 1.0e-6) const {
        if(!this->sameBins(o)) return false; //Check bins
        for(size_t d = 0; d < dim; ++d){ //Check limits
          //Handle near zero cases
          const double firstAbs = std::fabs(limits[d].first);
          const double secondAbs = std::fabs(limits[d].second);
          const double oFirstAbs = std::fabs(o.limits[d].first);
          const double oSecondAbs = std::fabs(o.limits[d].second);

          // + First limit check
          if(firstAbs < threshold || oFirstAbs < threshold){ // Near zero case
            if(firstAbs + oFirstAbs > 2.0*threshold) return false;
          }
          else{ //Non zero case
            const double relFirst  = limits[d].first/o.limits[d].first;
            if(std::fabs(relFirst - 1.0) > threshold) return false;
          }

          // + Second limit check
          if(secondAbs < threshold || oSecondAbs < threshold){ // Near zero case

            if(secondAbs + oSecondAbs > 2.0*threshold) return false;
            
          } else{ //Non zero case
            const double relSecond  = limits[d].second/o.limits[d].second;
            if(std::fabs(relSecond  - 1.0) > threshold) return false;
          }
        }
        return true;
      }

      // Print: same output your current printDims produces (low/top + width).
      void printDims(FILE* fout, const bool printOnlyEffective,
                     const std::array<std::string, dim>& headers) const {
        this->printBinning(fout, printOnlyEffective);
        fprintf(fout, "# Limits [low,top): \n");
        for(size_t d = 0; d < dim; ++d){
          if(printOnlyEffective && this->nBins[d] == 1) continue;
          fprintf(fout, "# %15.5E %15.5E  %s\n",
                  limits[d].first, limits[d].second,
                  headers[d].c_str());
        }
      }

      // Read the "# Limits [low,top):" block (low top header per line).
      // `headers` are filled by the caller from the trailing header text.
      int loadDims(std::istream& in,
                   std::array<std::string, dim>& headers,
                   unsigned long& readEffectiveDim) {

        unsigned long readDim;
        std::vector<unsigned long> readNBins;
        bool readEffective;
        int err = this->loadBinning(in, readDim, readNBins, readEffectiveDim, readEffective);
        if(err != errors::SUCCESS)
          return err;
        
        // ** Limits
        std::string line;
        bool found = false;
        while(std::getline(in, line)){
          if(line.find("# Limits [low,top):") == 0){
            found = true;
            break;
          }
        }

        if(!found){
          return errors::LIMITS_NOT_FOUND;
        }

        //Read limits
        size_t j = 0;
        for(size_t i = 0; i < readDim; ++i){
          if(readEffective && readNBins[i] == 1){
            //Skip line
            in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
          }
          else{
            //Read limits
            std::getline(in, line);

            //Create a buffer from this line
            std::stringstream ssline;
            ssline.str(line);

            //Skip comment character
            ssline.ignore(std::numeric_limits<std::streamsize>::max(), '#');

            //Read low limit
            ssline >> limits[j].first;
            if(!ssline){
              return errors::LIMITS_NOT_FOUND;
            }

            //Read top limit
            ssline >> limits[j].second;
            if(!ssline){
              return errors::LIMITS_NOT_FOUND;
            }

            //Read dimension header
            std::getline(ssline,headers[j]);

            if(headers[j].length() > 2){
              //Trim header right spaces
              headers[j] = trim(headers[j]);
            }
	    
            if(limits[j].first >= limits[j].second){
              return errors::INVALID_LIMITS;
            }
            ++j;
          }
        }

        //Calculate bin widths
        for(size_t i = 0; i < dim; ++i){
          binWidths[i] = (limits[i].second - limits[i].first)/static_cast<double>(this->nBins[i]);
        }

        return errors::SUCCESS;
      }
    };

    // ---------------------------------------------------------------------------
    //  VariableGrid: arbitrary (monotonically increasing) bin edges per dimension
    // ---------------------------------------------------------------------------
    template<size_t dim>
    class VariableGrid : public BinningBase<dim> {
    public:
      template<size_t newDim>
      using rebind = VariableGrid<newDim>;
      
      std::array<std::pair<double,double>, dim> limits{};
      std::array<std::vector<double>, dim> edges; // edges[d].size() == nBins[d]+1

      int init(const std::array<unsigned long, dim>& nBinsIn,
               const std::array<std::vector<double>, dim>& edgesIn) {
        for(size_t d = 0; d < dim; ++d){
          if(nBinsIn[d] == 0) return errors::INVALID_NUMBER_OF_BINS;
          if(edgesIn[d].size() != nBinsIn[d] + 1) return errors::DIMENSION_MISMATCH;
          for(size_t i = 1; i < edgesIn[d].size(); ++i)
            if(edgesIn[d][i] <= edgesIn[d][i-1]) return errors::INVALID_LIMITS;
        }
        this->nBins = nBinsIn;
        this->edges = edgesIn;
        this->computeIndexGeometry();
        for(size_t d = 0; d < dim; ++d){
          limits[d].first = edgesIn[d].front();
          limits[d].second = edgesIn[d].back();
        }
        return errors::SUCCESS;
      }

      std::array<unsigned long, dim>
      findBin(const std::array<double,dim>& pos) const {
        std::array<unsigned long, dim> idx;
        for(size_t d = 0; d < dim; ++d){
          auto it = std::upper_bound(edges[d].begin(), edges[d].end(), pos[d]);
          unsigned long i = static_cast<unsigned long>(it - edges[d].begin());
          idx[d] = (i == 0) ? 0 : std::min(i - 1, this->nBins[d] - 1);
        }
        return idx;
      }

      void bracket(const std::array<double,dim>& pos,
                   std::array<unsigned long,dim>& lo,
                   std::array<unsigned long,dim>& hi,
                   std::array<double,dim>& frac) const {
        for(size_t d = 0; d < dim; ++d){
          // Ensure the coordinate is within the domain [edges.front(), edges.back()]
          double x = pos[d];
          if(x < edges[d].front()){
            lo[d] = 0;
            hi[d] = 0;
            frac[d] = 0.0;
          }
          else if(x >= edges[d].back()){
            lo[d] = this->nBins[d]-1;
            hi[d] = this->nBins[d]-1;
            frac[d] = 1.0;
          }
          else{
            // Locate the bin: upper_bound returns first edge > x, so bin = it-1
            auto it = std::upper_bound(this->edges[d].begin(), this->edges[d].end(), x);
            unsigned long i = static_cast<unsigned long>(it - this->edges[d].begin());
            lo[d] = (i == 0) ? 0 : std::min(i - 1, this->nBins[d] - 1);

            // Upper neighbour
            hi[d] = std::min(lo[d] + 1, this->nBins[d] - 1);

            // Fraction is inside the bin: in [0,1)
            const double w = edges[d][lo[d] + 1] - edges[d][lo[d]];
            frac[d] = (x - edges[d][lo[d]]) / w;
          }
        }
      }      

      double binLower(size_t d, unsigned long i) const { return edges[d][i];     }
      double binUpper(size_t d, unsigned long i) const { return edges[d][i + 1]; }
      double lowEdge (size_t d) const { return edges[d].front(); }
      double highEdge(size_t d) const { return edges[d].back();  }
      double binWidth(size_t d, unsigned long i) const {
        return edges[d][i+1] - edges[d][i];
      }
      
      std::array<std::vector<double>, dim> getEdges(const bool last = true) const {
        if(last)
          return edges;
        else{
          std::array<std::vector<double>, dim> auxEdges = edges;
          for(std::vector<double>& v : auxEdges){
            v.pop_back();
          }
          return auxEdges;
        }
      }
      std::array<std::vector<double>, dim> getInitEdges() const { return edges; }
      template<size_t profDim>
      int getProfileInit(const std::array<size_t, profDim>& profDimsIndex,
                         const std::array<std::pair<unsigned long, unsigned long>, dim>& binLimits,
                         std::array<unsigned long, profDim>& profileBins,
                         std::array<std::vector<double>, profDim>& profEdges) const {

        int errBin = this->getProfileBins(profDimsIndex, binLimits, profileBins);
        if(errBin != errors::SUCCESS)
          return errBin;

        //Extract edges
        for(size_t i = 0; i < profDim; ++i){
          //Get dimension index in the original grid
          size_t profIndex = profDimsIndex[i];
          profEdges[i].resize(profileBins[i] + 1);
          for(unsigned long k = 0; k <= profileBins[i]; ++k){
            const unsigned long srcBin = binLimits[profIndex].first + k;
            profEdges[i][k] = this->binLower(profIndex, srcBin);
          }
        }
        
        return errors::SUCCESS;
      }
      
      const std::array<std::pair<double, double>, dim>& readLimits() const { return limits; }

      bool sameAs(const VariableGrid& o, const double threshold = 1.0e-5) const {
        if(!this->sameBins(o)) return false; //Check bins
        for(size_t d = 0; d < dim; ++d){ //Check dimensions
          if(edges[d].size() != o.edges[d].size()) return false;
          for(size_t e = 0; e < edges[d].size(); ++e){ //Check edges
            const double edgeAbs = std::fabs(edges[d][e]);
            const double oEdgeAbs = std::fabs(o.edges[d][e]);

            if(edgeAbs < threshold || oEdgeAbs < threshold){ //Near zero case
              if(edgeAbs + oEdgeAbs > 2.0*threshold) return false;
            }
            else{ //Non zero case
              const double rel = edges[d][e]/o.edges[d][e];
              if(std::fabs(rel - 1.0) > threshold) return false;
            }
          }
        }
        return true;
      }

      // Variable grid prints its edge list instead of low/top + width.
      void printDims(FILE* fout, const bool printOnlyEffective,
                     const std::array<std::string, dim>& headers) const {
        this->printBinning(fout, printOnlyEffective);
        fprintf(fout, "# Bin edges:\n");
        for(size_t d = 0; d < dim; ++d){
          if(printOnlyEffective && this->nBins[d] == 1) continue;
          fprintf(fout, "#");
          for(double e : edges[d]) fprintf(fout, " %15.5E", e);
          fprintf(fout, "  %s\n", headers[d].c_str());
        }
      }

      // Reads nBins[d]+1 doubles per line.
      int loadDims(std::istream& in,
                   std::array<std::string, dim>& headers,
                   unsigned long& readEffectiveDim)  {

        unsigned long readDim;
        std::vector<unsigned long> readNBins;
        bool readEffective;
        int err = this->loadBinning(in, readDim, readNBins, readEffectiveDim, readEffective);
        if(err != errors::SUCCESS)
          return err;
        
        // ** Limits
        std::string line;
        bool found = false;
        while(std::getline(in, line)){
          if(line.find("# Bin edges:") == 0){
            found = true;
            break;
          }
        }

        if(!found){
          return errors::LIMITS_NOT_FOUND;
        }

        //Set default edges
        for(size_t d = 0; d < dim; d++){
          edges[d].resize(2);
          edges[d][0] = -1.0e35;
          edges[d][1] = +1.0e35;
        }

        //Read limits
        size_t j = 0;
        for(size_t i = 0; i < readDim; ++i){
          if(readEffective && readNBins[i] == 1){
            //Skip line
            in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
          }
          else{
            //Read limits
            std::getline(in, line);

            //Create a buffer from this line
            std::stringstream ssline;
            ssline.str(line);

            //Skip comment character
            ssline.ignore(std::numeric_limits<std::streamsize>::max(), '#');

            //Read edges
            edges[j].resize(this->nBins[j]+1);
            for(double& edge : edges[j]){
              ssline >> edge;
              if(!ssline){
                return errors::LIMITS_NOT_FOUND;
              }
            }

            //Read dimension header
            std::getline(ssline,headers[j]);

            if(headers[j].length() > 2){
              //Trim header right spaces
              headers[j] = trim(headers[j]);
            }
	    
            //Check edges
            for(size_t k = 1; k < edges[j].size(); ++k){
              if(edges[j][k] <= edges[j][k-1]) return errors::INVALID_LIMITS;
            }

            ++j;
          }
        }

        //Set limits
        for(size_t d = 0; d < dim; ++d){
          limits[d].first = edges[d].front();
          limits[d].second = edges[d].back();
        }

        return errors::SUCCESS;
      }
    };

  } // namespace grids
} // namespace penred
#endif
