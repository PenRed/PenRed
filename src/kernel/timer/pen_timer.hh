
//
//
//    Copyright (C) 2019-2020 Universitat de València - UV
//    Copyright (C) 2019-2020 Universitat Politècnica de València - UPV
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

 
#ifndef __PEN_TIMER__
#define __PEN_TIMER__

#include <ctime>
#include <cstring>
#include <cstdio>
#include <chrono>

class pen_timer{

private:
  time_t time_beginning;
public:


  //  ****  A call to TIMEO initialises the clock.
  inline void time0(){time_beginning = time(NULL);} 

  double timer();
  pen_timer(){time0();}
  
};

class pen_stopWatch{
private:
  std::chrono::steady_clock::time_point stop;
  std::chrono::milliseconds interval;
public:

  pen_stopWatch() : interval(0){
    stop = std::chrono::steady_clock::now();
  }

  template<class T> pen_stopWatch(T elapsed, const typename std::enable_if<std::is_integral<T>::value,bool>::type /*dummy*/ = true){
    std::chrono::milliseconds::rep toWait;
    if(std::numeric_limits<std::chrono::milliseconds::rep>::max() < elapsed)
      toWait = std::numeric_limits<std::chrono::milliseconds::rep>::max();
    else
      toWait = static_cast<std::chrono::milliseconds::rep>(elapsed);
    
    interval = std::chrono::milliseconds(toWait);
    stop = std::chrono::steady_clock::now();
  }
  
  template<class T> void duration(const T elapsed){
    interval = std::chrono::milliseconds(elapsed);
  }

  inline void start(){
    stop = std::chrono::steady_clock::now();
    stop += interval;
  }

  inline bool check(const std::chrono::steady_clock::time_point& time) const {
    if(time >= stop)
      return true;
    return false;
  }
  
  inline bool check() const {
    if(std::chrono::steady_clock::now() >= stop)
      return true;
    return false;
  }
};
  
double CPUtime();
void PDATET(char* DATE23);

#endif
