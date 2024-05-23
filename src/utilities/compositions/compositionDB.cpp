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
//        vicent.gimenez.alventosa@gmail.com (Vicent Giménez Alventosa)
//    
//

#include "composDB.hh"

int main (int argc, const char** argv){

  printf("usage: %s [composition-DB] [material-name]\n", argv[0]);

  if(argc == 1){
    printf("Registered data bases:\n");
    for(const char* s : penred::dataBases::compositions::dataBases){
      printf(" + %s\n", s);
    }
    return 0;
  }

  //Get composition DB pointer
  auto db = penred::dataBases::compositions::getDB(argv[1]);

  if(db == nullptr){
    printf("Unknown composition database '%s'\n", argv[1]);
    printf("Registered data bases:\n");
    for(const char* s : penred::dataBases::compositions::dataBases){
      printf(" + %s\n", s);
    }
    return 1;    
  }

  //If no material has been selected, print all availables
  if(argc == 2){
    printf("Registered compositions in database '%s':\n", argv[1]);
    for(const std::string& s : db->matList()){
      printf("  + %s\n", s.c_str());
    }
    return 0;
  }

  //Print the specified materials
  for(int imat = 2; imat < argc; ++imat){

    const std::vector<penred::massFraction> mfs = db->getElements(argv[imat]);
    const double dens = db->getDensity(argv[imat]);
    if(mfs.size() == 0 || dens <= 0.0){
      printf("\n  - Unknown material '%s'.\n", argv[imat]);
    }
    else{
      printf("\n  + Material '%s' (%15.5E g/cm**3):\n",
	     argv[imat], dens);

      for(size_t j = 0; j < mfs.size(); ++j){
	printf("    - %3u  %f\n", mfs[j].Z, mfs[j].fraction);
      }
    }
  }

  return 0;
}
