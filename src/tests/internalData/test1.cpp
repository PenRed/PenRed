 
//
//
//    Copyright (C) 2019 Universitat de València - UV
//    Copyright (C) 2019 Universitat Politècnica de València - UPV
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

 
 
#include "pen_data.hh"

int main(){

    const char c = 'u';
    const int i = 12345;
    const double d = 1.53e6;
    const bool b = true;
    
    pen_parserSection root;
    pen_parserArray array;

    pen_parserData data;
    pen_parserData data1(c);
    pen_parserData data2(i);
    pen_parserData data3(d);
    pen_parserData data4(b);
    
    array.append(data1);
    array.append(data2);
    array.append(data3);
    
    int err;
    err = root.set("data4",data4);
    if(err != INTDATA_SUCCESS){
        printf("set1:Error code %d\n",err);
    }
        
    err = root.set("array1",array);
    if(err != INTDATA_SUCCESS){
        printf("set2:Error code %d\n",err);
    }

    err = root.set("futureSection",546);
    if(err != INTDATA_SUCCESS){
        printf("set3:Error code %d\n",err);
    }

    err = root.set("sec0/element32",98.2e2);
    if(err != INTDATA_SUCCESS){
        printf("set4:Error code %d\n",err);
    }
    
    err = root.set("sec0/sec1/en","hi");
    if(err != INTDATA_SUCCESS){
        printf("set5:Error code %d\n",err);
    }
    
    std::string auxstr("sec0/sec1/es");
    err = root.set(auxstr,"hola");
    if(err != INTDATA_SUCCESS){
        printf("set6:Error code %d\n",err);
    }
    
    err = root.set("sec0/sec1/¿?",1234);
    if(err != INTDATA_SUCCESS){
        printf("set7:Error code %d\n",err);
    }
    
    err = root.set("sec0/sec2/array_cp",array);
    if(err != INTDATA_SUCCESS){
        printf("set8:Error code %d\n",err);
    }
    err = root.set("sec0/sec2/array_cp2",array);
    if(err != INTDATA_SUCCESS){
        printf("set9:Error code %d\n",err);
    }

    err = root.set("sec6/scalar",'v');
    if(err != INTDATA_SUCCESS){
        printf("set10:Error code %d\n",err);
    }
    
    
    //Stringify root section
    std::string rootString;
    root.stringify(rootString);
    printf("%s\n",rootString.c_str());
    
    //Get section /sec0/sec1
    pen_parserSection sec0sec1;
    err = root.readSubsection("/sec0/sec1",sec0sec1);
    if(err != INTDATA_SUCCESS){
        printf("readSubsection:Error code %d\n",err);
    }

    printf("\n\n");
    printf("---------------------------------------------\n");
    
    //Get ls to /sec0/sec1
    std::vector<std::string> vectSec0;
    err = root.ls("/sec0/",vectSec0);
    if(err != INTDATA_SUCCESS){
        printf("ls:Error code %d\n",err);
    }

    for(unsigned j = 0; j < vectSec0.size(); j++){
      printf("%s\n",vectSec0[j].c_str());
    }

    printf("\n\n");
    vectSec0.clear();
    err = root.ls("/sec0",vectSec0);
    if(err != INTDATA_SUCCESS){
        printf("ls:Error code %d\n",err);
    }

    for(unsigned j = 0; j < vectSec0.size(); j++){
      printf("%s\n",vectSec0[j].c_str());
    }

    printf("\n\n");
    vectSec0.clear();
    err = root.ls("/sec0/element32",vectSec0);
    if(err != INTDATA_SUCCESS){
        printf("ls:Error code %d\n",err);
    }

    for(unsigned j = 0; j < vectSec0.size(); j++){
      printf("%s\n",vectSec0[j].c_str());
    }        
    //Print subsection
    printf("---------------------------------------------\n");
    std::string sec0sec1String;
    sec0sec1.stringify(sec0sec1String);
    printf("%s\n",sec0sec1String.c_str());
    
    return 0;
}
