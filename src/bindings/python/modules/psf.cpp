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

	  //Check parameters
	  if(emin >= emax){
	    throw py::value_error("Error:'emin' must be greater than 'emax'");
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
    emin (double): Lower bound of energy range in eV (must be >= 0).
emax (double): Upper bound of energy range in eV (must be > emin).
nbins (unsigned): Number of linear-spaced energy bins between emin and emax.

Returns:
    tuple: A tuple containing four tuples in the following order: energy bin edges, electron, gamma and positron spectrum.
		   For example, for a emin = E1, emax = En, and n bins, returns:
			((E1, E2, ..., En),(e- (E1), ..., e- (En)), (gamma (E1), ..., gamma (En)), (e+ (E1), ..., e+ (En)))

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
