//
//
//    Copyright (C) 2024 Universitat de València - UV
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
//        vicent.gimenez.alventosa@gmail.com (Vicent Giménez Alventosa)
//        sanolgi@upvnet.upv.es (Sandra Oliver Gil)
//    
//

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "functions.hh"
#include "pen_phaseSpaceFile.hh"


namespace py = pybind11;


PYBIND11_MODULE(psf,m){

  m.doc() = "penred phase space file module";
  
  m.def("psfSpectre",
        [](const std::string& filename,
           const double emin,
           const double emax,
           const unsigned nbins)->py::tuple{

	  //Check parameters
	  if(emin < 0){
	    throw py::value_error("Error:'emin' must be >= 0");
	  }

	  if(emin >= emax){
	    throw py::value_error("Error:'emin' must be greater than 'emax'");
	  }

	  if(nbins == 0){
	    throw py::value_error("Error: At least one energy been is necessary");
	  }
      
	  //Calculate bin width
	  double de = (emax - emin) / (double)nbins;
	  double ide = 1.0/de;

	  //Create a phase space file
	  pen_psfreader psf;


	  //Open specified input file
	  FILE* fin = nullptr;
	  fin = fopen(filename.c_str(),"rb");
	  if(fin == nullptr){
	    throw py::value_error("Error: unable to open file: " + filename);
	  }

	  //Create tally variables
	  std::array<std::vector<double>, constants::nParTypes> spectre;
	  for(std::vector<double>& s : spectre){
	    s.resize(nbins, 0.0);
	  }	  

	  //Read input file until the end
	  unsigned nchunks = 0;
	  long long unsigned nhists = 0;
	  while(psf.read(fin,1) == PEN_PSF_SUCCESS){
	    nchunks++;
	    //Iterate over read states
	    pen_particleState state;
	    unsigned long dhist;
	    unsigned kpar;
	    while(psf.get(dhist,kpar,state) > 0){

	      //Calculate bin index
	      int ibin = (state.E-emin)*ide;

	      nhists += dhist;
	      //Check bin index
	      if(ibin >= 0 && ibin < (int)nbins){
		spectre[kpar][ibin] += state.WGHT;
	      }
	    }
	  }

	  //Close files
	  fclose(fin);
		

	  //Create the resulting python tuple
	  py::tuple results(constants::nParTypes+1);

	  std::vector<double> energies(nbins);
	  for(unsigned j = 0; j < nbins; j++)
	    {
	      energies[j] = double(j)*de+emin;
	    }

	  for(unsigned i = 0; i < constants::nParTypes; i++)
	    {
	      std::vector<double> spectPSF(nbins);
	      for(unsigned j = 0; j < nbins; j++)
		{
		  spectPSF[j] = spectre[i][j];
		}
	      results[i+1] = py::tuple(py::cast(spectPSF));
	    }

		
	  results[0] = py::tuple(py::cast(energies));
	  return results;
	},    
	py::arg("filename"),
	py::arg("emin"),
	py::arg("emax"),
	py::arg("nbins"),
	R"doc(
Extracts the energy spectrum distribution of the particles from a Phase Space File.

Args:
    filename (str): Path to the Phase Space File.
    emin (float): Lower bound of energy range in eV (must be >= 0).
    emax (float): Upper bound of energy range in eV (must be > emin).
    nbins (unsigned): Number of linear-spaced energy bins between emin and emax.

Returns:
    tuple: A tuple containing four tuples in the following order: energy bin edges, electron, gamma and positron spectrum.

Example:
    .. code-block:: python

		# Read the Phase space File and store the results
			results = pyPenred.psf.psfSpectre('psf.dat', 0, 7e5, 200)

)doc");

  m.def("getDistribution",
        [](const std::string& filename,
           const double emin,
           const double emax,
           const unsigned ebins,
           const double xmin,
           const double xmax,
           const unsigned xbins,
           const double ymin,
           const double ymax,
           const unsigned ybins,
           const double zmin,
           const double zmax,
           const unsigned zbins,
           const unsigned cosBins,
           const bool toFile,
           const bool extractInfo,
           const bool onlyEffective)->py::tuple{

          //Check parameters
          if(emin < 0){
            throw py::value_error("Error:'emin' must be >= 0");
          }

          if(emin >= emax){
            throw py::value_error("Error:'emin' must be greater than 'emax'");
          }
          
          if(xmin >= xmax){
            throw py::value_error("Error:'xmin' must be greater than 'xmax'");
          }

          if(ymin >= ymax){
            throw py::value_error("Error:'ymin' must be greater than 'ymax'");
          }

          if(zmin >= zmax){
            throw py::value_error("Error:'zmin' must be greater than 'zmax'");
          }

          if(cosBins == 0){
            throw py::value_error("Error: At least one cosinus bin is necessary");
          }
          if(ebins == 0){
            throw py::value_error("Error: At least one energy bin is necessary");
          }
          if(xbins == 0 || ybins == 0 || zbins == 0){
            throw py::value_error("Error: At least one spatial been in each axis is necessary");
          }

          
                    
          //Create a phase space file
          pen_psfreader psf;

          //Open specified input file
          FILE* fin = nullptr;
          fin = fopen(filename.c_str(),"rb");
          if(fin == nullptr){
            throw py::value_error("Error: unable to open file: " + filename);
          }

          //Create the tallies
          penred::measurements::measurement<double, 4> energyDistrib;
          energyDistrib.description = "PSF energetic spatial distribution";
  
          energyDistrib.setDimHeader(0, "E (eV)");
          energyDistrib.setDimHeader(1, "x (cm)");
          energyDistrib.setDimHeader(2, "y (cm)");
          energyDistrib.setDimHeader(3, "z (cm)");
          energyDistrib.setValueHeader("Prob(1/hist)");

          energyDistrib.initFromLists({ebins, xbins, ybins, zbins},
                                      {penred::measurements::limitsType(emin, emax),
                                       penred::measurements::limitsType(xmin, xmax),
                                       penred::measurements::limitsType(ymin, ymax),
                                       penred::measurements::limitsType(zmin, zmax)});          
          

          penred::measurements::measurement<double, 4> cosDistrib;
          cosDistrib.description = "Direction cosinus spatial distribution. It is calculated respect Z axis";
  
          cosDistrib.setDimHeader(0, "cos");
          cosDistrib.setDimHeader(1, "x (cm)");
          cosDistrib.setDimHeader(2, "y (cm)");
          cosDistrib.setDimHeader(3, "z (cm)");
          cosDistrib.setValueHeader("Prob(1/hist)");
          
          cosDistrib.initFromLists({cosBins, xbins, ybins, zbins},
                                   {penred::measurements::limitsType(-1.0, 1.000001),
                                    penred::measurements::limitsType(xmin, xmax),
                                    penred::measurements::limitsType(ymin, ymax),
                                    penred::measurements::limitsType(zmin, zmax)});          
          
          //Read input file until the end
          unsigned nchunks = 0;
          long long unsigned nhist = 0;
          while(psf.read(fin,1) == PEN_PSF_SUCCESS){
            nchunks++;
            //Iterate over read states
            pen_particleState state;
            unsigned long dhist;
            unsigned kpar;
            while(psf.get(dhist,kpar,state) > 0){
              nhist += dhist;
              energyDistrib.add({state.E, state.X, state.Y, state.Z}, state.WGHT, nhist);
              cosDistrib.add({state.W, state.X, state.Y, state.Z}, state.WGHT, nhist);
            }
          }

          //Close files
          fclose(fin);

          //Create the resulting python tuples
          penred::measurements::results<double, 4> resultE;
          energyDistrib.results(nhist, resultE);

          penred::measurements::results<double, 4> resultCos;
          cosDistrib.results(nhist, resultCos);
          
          //Check if the distribution shoud be saved
          if(toFile){
            FILE* fout = nullptr;
            fout = fopen("psf-eDistribution.dat", "w");
            if(fout == nullptr){
              throw std::runtime_error("Error: Unable to open file 'psf-eDistribution.dat'");
            }
            resultE.print(fout, 2, true, true, onlyEffective);
            fclose(fout);
            fout = nullptr;

            if(cosBins > 1){
              fout = fopen("psf-cosDistribution.dat", "w");
              if(fout == nullptr){
                throw std::runtime_error("Error: Unable to open file 'psf-cosDistribution.dat'");
              }
              resultCos.print(fout, 2, true, true, onlyEffective);
              fclose(fout);
            }
          }
          
          py::tuple resTupleE = result2numpy(resultE, extractInfo, onlyEffective);
          py::tuple resTupleCos = result2numpy(resultCos, extractInfo, onlyEffective);

          if(cosBins > 1)
            return py::make_tuple(resTupleE,resTupleCos);
          else
            return resTupleE;
        },
        py::arg("filename"),
        py::arg("emin") = 0.0,
        py::arg("emax") = 1.0e9,
        py::arg("ebins") = 1,
        py::arg("xmin") = -1.0e30,
        py::arg("xmax") =  1.0e30,
        py::arg("xbins") = 1,
        py::arg("ymin") = -1.0e30,
        py::arg("ymax") =  1.0e30,
        py::arg("ybins") = 1,
        py::arg("zmin") = -1.0e30,
        py::arg("zmax") =  1.0e30,
        py::arg("zbins") = 1,
        py::arg("cos_bins") = 1,
        py::arg("save_files") = false,
        py::arg("extract_info") = true,        
        py::arg("only_effective") = true,        
        R"doc(
Extracts the energy spectrum distribution of the particles from a Phase Space File.

Args:
    filename (str): Path to the Phase Space File.
    emin (float): Lower energy bound in eV (must be >= 0).
    emax (float): Upper energy bound in eV (must be > emin).
    ebins (unsigned): Number of linear-spaced energy bins between emin and emax.
    xmin (float): Lower spatial bound along the X-axis.
    xmax (float): Upper spatial bound along the X-axis.
    xbins (unsigned): Number of linear-spaced spatial bins along the X-axis.
    ymin (float): Lower spatial bound along the Y-axis.
    ymax (float): Upper spatial bound along the Y-axis.
    ybins (unsigned): Number of linear-spaced spatial bins along the Y-axis.
    zmin (float): Lower spatial bound along the Z-axis.
    zmax (float): Upper spatial bound along the Z-axis.
    zbins (unsigned): Number of linear-spaced spatial bins along the Z-axis.
    cos_bins (unsigned): Number of linear-spaced cosinus bins between -1 and 1.
    save_files (bool): If enabled, the distribution is saved in a file named 'psf-distribution.dat'.
    extract-info (bool): If enabled, the limits and dimensions information will be returned along with results values.
    only_effective (bool): If enabled, dimensions with a single bin are ignored, reducing the result's dimension.

Returns:
    tuple: A tuple containing four tuples in the following order: energy bin edges, electron, gamma and positron spectrum.

Example:
    .. code-block:: python

		# Read the Phase space File and store the results
			results = pyPenred.psf.psfSpectre('psf.dat', 0, 7e5, 200)

)doc");

m.def("psf2ascii",
      [](const std::string& filenameIn,
	 const std::string& filenameOut)->py::tuple{
		   
	if(filenameIn == filenameOut)
	  throw py::value_error("Error: output ASCII filename (" + filenameOut + ") must be different from input binary filename (" + filenameIn + ")");

	//Create a phase space file
	pen_psfreader psf;	
		

	//Open specified input file
	FILE* fin = nullptr;
	fin = fopen(filenameIn.c_str(),"rb");
	if(fin == nullptr){
	  throw py::value_error("Error: unable to open input file: " + filenameIn);
	}
		
	//Open output file
	FILE* fout = nullptr;
	fout = fopen(filenameOut.c_str(),"w");
	if(fout == nullptr){
	  throw py::value_error("Error: unable to open output file: " + filenameOut);
	}
		
	//Print header
	fprintf(fout,"# DHIST  KPAR %s\n",baseStateHeaderNoGeo());

	//Read input file until the end
	unsigned nchunks = 0;
	long long unsigned nhists = 0;
	while(psf.read(fin,1) == PEN_PSF_SUCCESS){
	  nchunks++;
	  //Iterate over read states
	  pen_particleState state;
	  unsigned long dhist;
	  unsigned kpar;
	  while(psf.get(dhist,kpar,state) > 0){
	    nhists += dhist;
	    //Print base state without body and material information
	    fprintf(fout,"  %5lu  %4u %s\n",dhist,kpar,state.stringifyBaseNoGeo().c_str());
	  }
	}

	//Close files
	fclose(fin);
	fclose(fout);  

	//Create the resulting python tuple
	py::tuple results(2);

	results[0]= py::cast(nhists);
	results[1]= py::cast(nchunks);

	return results;
      },

		py::arg("filenameIn"),
		py::arg("filenameOut"),
		R"doc(
	Converts the binary Phase Space File into an ASCII file format.

	Args:
		filenameIn (str): Path to the Phase Space File.
		filenameOut (str): Path to the converted Phase Space File into an ASCII format.

	Returns:
		tuple: A tuple containing the number of histories in the Phase Space File and the number of particle chuncks.

	Example:
		.. code-block:: python

			# Read the Phase space File and store the results
				results = pyPenred.psf.psf2ascii('psf.dat', 'psfASCII.dat')

)doc");
}
