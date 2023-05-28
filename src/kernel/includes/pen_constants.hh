
//
//
//    Copyright (C) 2019-2023 Universitat de València - UV
//    Copyright (C) 2019-2023 Universitat Politècnica de València - UPV
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


#ifndef __PEN_CONSTANTS_
#define __PEN_CONSTANTS_

#include "../particles/includes/pen_particles_ID.hh"

namespace constants{

const double A0B = 5.2917721092E-9;  // Bohr radius (cm)
const double HREV = 27.21138505;  // Hartree energy (eV)
const double AVOG = 6.02214129E23;  // Avogadro's number
const double SL = 137.035999074;  // Speed of light (1/alpha)
const double PI = 3.1415926535897932; // PI
const double TWOPI = 2.0*PI;
const double REV = 5.10998928E5;  // Electron rest energy (eV)
const double ELRAD = 2.8179403267E-13;  // Class. electron radius (cm)
const unsigned int nParTypes = pen_KPAR::ALWAYS_AT_END; //Number of particle types
const unsigned int NO = 512; //****  E/P inelastic collisions maximum number of oscillators.
const unsigned int MAXMAT = 100; //Maximum number of materials
const unsigned int MAXINTERACTIONS = 8; //Maximum number of interactions per particle type
const unsigned int NEGP = 200; //Energy grid points
const unsigned int NBW = 32; //  ****  Bremsstrahlung emission.
const unsigned int NBE = 57;

const unsigned int NOCO = 512; //  ****  Compton scattering. number of oscillators.

//  ****  Rayleigh scattering.
const unsigned int NQ = 250;
const unsigned int NEX = 1024;
 
const unsigned int NTP = 15000; //  ****  Photoelectric cross sections.

const unsigned long int NRX = 60000; // RELAX

const unsigned int NMS = 5000; //Secondary stack size

const double MINEABSE = 50.0; // **** Minima energia d'absorcio en electrons (eV)
const double MINEABSPh = 50.0; // **** Minima energia d'absorcio en fotons (eV)
const double MINEABSPo = 50.0; // **** Minima energia d'absorcio en positrons (eV)
const double MINEGRID = 50.0; // **** Minima energia del "grid" (eV). Deuria ser la minima de les energies d'absorcio.
const double MINDIMGRID = 10.0; // **** Minima distancia entre el valor maxim i minim del "grid" (eV)

const double HBAR = 6.58211928E-16;

// Electrons hard bremsstrahlung emission
const double WB[constants::NBW] = {1.0E-12, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.925, 0.95, 0.97, 0.99, 0.995, 0.999, 0.9995, 0.9999, 0.99995, 0.99999, 1.0};
  
}

namespace pen_geoconst{
  const unsigned int NB  = 5000;  //Maximum number of bodies in geometry
}

#endif
