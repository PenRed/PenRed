
//
//
//    Copyright (C) 2019-2020 Universitat de València - UV
//    Copyright (C) 2019-2020 Universitat Politècnica de València - UPV
//    Copyright (C) 2024-2025 Vicent Giménez Alventosa
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
  std::chrono::steady_clock::time_point _start;
  std::chrono::steady_clock::time_point stop;
  std::chrono::milliseconds interval;
  
public:

  static constexpr const std::chrono::milliseconds MAX_DURATION_MS{157'788'000'000}; //Allow "only" 5 years as time intervals
  
  pen_stopWatch() : interval(0){
    _start = std::chrono::steady_clock::now();
    stop = _start;
  }

  template<class T> pen_stopWatch(T elapsed, const typename std::enable_if<std::is_integral<T>::value,bool>::type /*dummy*/ = true){

    if(elapsed > 0){    
      if(static_cast<double>(MAX_DURATION_MS.count()) <= static_cast<double>(elapsed))
	interval = MAX_DURATION_MS;
      else
	interval = std::chrono::milliseconds(static_cast<std::chrono::milliseconds::rep>(elapsed));
    }else{
      interval = std::chrono::milliseconds(0);
    }
    _start = std::chrono::steady_clock::now();
    stop = _start;
  }
  
  template<class T> void duration(const T elapsed){

    if(elapsed > 0){    
      if(static_cast<double>(MAX_DURATION_MS.count()) <= static_cast<double>(elapsed))
	interval = MAX_DURATION_MS;
      else
	interval = std::chrono::milliseconds(static_cast<std::chrono::milliseconds::rep>(elapsed));
    }else{
      interval = std::chrono::milliseconds(0);
    }

    if(_start < stop){
      //Update stop
      stop = _start + interval;
    }
  }

  inline void start(){
    _start = std::chrono::steady_clock::now();
    stop = _start + interval;
  }

  inline bool check(const std::chrono::steady_clock::time_point& time) const {
    return time >= stop;
  }
  
  inline bool check() const {
    return std::chrono::steady_clock::now() >= stop;
  }
};
  
double CPUtime();
void PDATET(char* DATE23);

#endif
