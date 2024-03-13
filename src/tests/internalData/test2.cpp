 
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
 
    char cRead = 'o';
    int iRead = 0;
    double dRead = 0.0;
    bool bRead = 0;
    
    pen_parserData data;
    pen_parserData data1(c);
    pen_parserData data2(i);
    pen_parserData data3(d);
    pen_parserData data4(b);
    
    
    //Check read
    printf("READ CHECK ----------------------------\n");
    
    int err;
    err = data1.read(iRead);
    if(err != INTDATA_SUCCESS)
    {
        printf("set1: Error code %d\n", err);
    }
    else
    {
      printf("set1: cRead %c == %s\n", cRead, data1.stringify().c_str());
    }
    
    err = data2.read(iRead);
    if(err != INTDATA_SUCCESS)
    {
        printf("set2: Error code %d\n", err);
    }
    else
    {
      printf("set2: iRead %i == %s\n", iRead, data2.stringify().c_str());
    }
    
    err = data3.read(dRead);
    if(err != INTDATA_SUCCESS)
    {
        printf("set3: Error code %d\n", err);
    }
    else
    {
        printf("set3: dRead %lf == %s\n", dRead, data3.stringify().c_str());
    }
    

    err = data4.read(bRead);
    if(err != INTDATA_SUCCESS)
    {
        printf("set4: Error code %d\n", err);
    }
    else
    {
        printf("set4: bRead %s == %s\n", (bRead) ? "true" : "false", data4.stringify().c_str());
    }
    
    printf("---------------------------------------\n");
    
    
    
    
    //-----------------------------------------------------//
    
    
    //Check readTag
    err = data1.readTag();
    if(err != CHAR)
    {
        printf("set5: Error code %d\n", err);
    }
    
    err = data2.readTag();
    if(err != INT)
    {
        printf("set6: Error code %d\n", err);
    }
    
    err = data3.readTag();
    if(err != DOUBLE)
    {
        printf("set7: Error code %d\n", err);
    }
    
    err = data4.readTag();
    if(err != BOOL)
    {
        printf("set8: Error code %d\n", err);
    }
   
   
 
    //printf(" data 1 tag = %d\n data 2 tag = %d\n data 3 tag = %d\n data 4 tag = %d\n", data1.readTag(), data2.readTag(), data3.readTag(), data4.readTag());
    
   
    
}
