
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


#ifndef __SUBCLASS_INSTANTIATOR__
#define __SUBCLASS_INSTANTIATOR__

#include <cstring>
#include <mutex>
#include <functional>
#include <vector>
#include <map>

template <class mother, class sub>
void instanceInheritance(mother*& pmoth){
  sub* psub = new sub;
  pmoth = nullptr;
  pmoth = dynamic_cast<mother*>(psub);
  if(pmoth == nullptr){
    delete psub;
  }
}

template <class motherClass>
class instantiator
{
  typedef std::map<std::string,std::function<void (motherClass*&)>> mapType;
private:

  std::mutex lock;
  mapType creatorsMap;
  
public:

  instantiator()
  {}

  inline unsigned registred(){
    return creatorsMap.size();
  }

  template <class subclass>
  int addSubType(const char* typeID){
    return addSubType<subclass>(std::string(typeID));
  }
  template <class subclass>
  int addSubType(const std::string& typeID){
    
    //Lock the mutex
    std::lock_guard<std::mutex> guard(lock);
    
    //Check if this ID has been already registered
    if(creatorsMap.find(typeID) != creatorsMap.end())
      return 1;
    
    //Check if is possible to instance mother class from subclass
    motherClass* pclass;
    instanceInheritance<motherClass,subclass>(pclass);
    
    if(pclass != nullptr){

      //Free memory
      delete pclass;
      
      //Get function 
      std::function<void (motherClass*&)> func = instanceInheritance<motherClass,subclass>;
    
      //register new subclass
      creatorsMap[typeID] = func;

      //Return success
      return 0;
    }
    //Return fail
    return -1;
  }

  motherClass* createInstance(const char* typeID) {
    return createInstance(std::string(typeID));
  }
  motherClass* createInstance(std::string typeID) {

    //Lock the mutex
    std::lock_guard<std::mutex> guard(lock);

    //Search specified ID
    if(creatorsMap.find(typeID) != creatorsMap.end()){
      motherClass* pmoth = 0;
      creatorsMap[typeID](pmoth);
      return pmoth;
    }
    
    return nullptr;
  }

  std::string typesList(){

    //Lock the mutex
    std::lock_guard<std::mutex> guard(lock);

    std::string aux;
    typename mapType::const_iterator it;
    for(it = creatorsMap.begin(); it != creatorsMap.end(); it++){
      aux += it->first + "\n";
    }
    return aux;
  }

  std::string typesList(std::vector<std::string>& list){

    //Lock the mutex
    std::lock_guard<std::mutex> guard(lock);

    std::string aux;
    typename mapType::const_iterator it;
    for(it = creatorsMap.begin(); it != creatorsMap.end(); it++){
      aux += it->first + "\n";
      list.push_back(it->first);
    }
    return aux;
  }  
  
};


#endif
