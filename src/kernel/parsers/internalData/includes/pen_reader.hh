
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
//        vicent.gimenez.alventosa@gmail.com
//        vicente.gimenez@uv.es
//    
//

 
#ifndef __PENRED_READER_CLASSES__
#define __PENRED_READER_CLASSES__

#include "pen_parser.hh"

//Template wrapper to store specific read formats for each configurable object
template<class configurableObj>
struct pen_format{
  //By default, no format is defined
  static constexpr const char* format = nullptr;  
};

class pen_readerSection;

struct pen_readerRequired{

  static const std::string typeKey;
  static const std::string valueKey;
  
  enum requiredTypes{
    REQUIRED,
    REQUIRED_IF,
    REQUIRED_IF_EXIST,
    OPTIONAL,
    OPTIONAL_IF,
    OPTIONAL_IF_EXIST,
  };

  requiredTypes reqType;
  std::string reqPath;
  
  pen_readerRequired() : reqType(REQUIRED){}

  //Parse function
  int parse(const pen_parserSection& s, const std::string& prefixIn = "");
  
  bool required(const pen_parserSection& s,
		const std::string& actualSectionPath = "") const;

  inline std::string stringify() const {
    std::string result;

    switch(reqType){
    case REQUIRED:     return "value is required";
    case REQUIRED_IF:{
      result = "value is required if '" + reqPath + "' is true" ;
      return result;
    }
    case REQUIRED_IF_EXIST:{
      result = "value is required if '" + reqPath + "' exists" ;
      return result;
    }
    case OPTIONAL:     return "value is optional";
    case OPTIONAL_IF:{
      result = "value is optional if '" + reqPath + "' is true" ;
      return result;
    }
    case OPTIONAL_IF_EXIST:{
      result = "value is optional if '" + reqPath + "' exists" ;
      return result;
    }
    }
    
    return "value is required";
  }
  inline std::string stringify(const pen_parserSection& s) const {

    std::string result = stringify();

    if(required(s))
      result.append(" (Required) ");
    else
      result.append(" (Optional) ");

    return result;
  }

  inline void clear(){
    reqType = REQUIRED;
    reqPath.clear();
  }
};

struct pen_readerCondition{

  static const std::string typeKey;
  static const std::string descriptionKey;
  static const std::string valueKey;
  
  enum conditionTypes{
    POSITIVE,
    NEGATIVE,
    LESSER,
    LESSER_EQUAL,
    GREATER,
    GREATER_EQUAL,
    NOT_EQUAL,
    EQUAL,
    NO_CONDITION,
  };

  enum valueTypes{
    NONE,
    INT,
    DOUBLE,
    PATH,
  };
  
  conditionTypes condType;
  valueTypes valType;

  std::string description;
  std::string condPath;
  int condInt;
  double condDouble;

private:

  template<class Tread, class Tcond>
  bool checkT(const Tread readValue, const Tcond checkValue) const {
    
    //Convert read type to condition type
    Tcond value = static_cast<Tcond>(readValue);

    switch(condType){
    case LESSER:        return value < checkValue;
    case LESSER_EQUAL:  return value <= checkValue;
    case GREATER:       return value > checkValue;
    case GREATER_EQUAL: return value >= checkValue;
    case NOT_EQUAL:     return value != checkValue;
    case EQUAL:         return value == checkValue;
    case NO_CONDITION:  return true;
    default: return true;
    }
    
    //No condition
    return true;
  }
  
public:
  
  pen_readerCondition() : condType(NO_CONDITION),
			  valType(NONE)
  {}

  //Parse function
  int parse(const pen_parserSection& s,
	    std::string& errorString,
	    const std::string& prefixIn = "");
  
  template<class T>
  typename std::enable_if<std::is_arithmetic<T>::value, bool>::type
  check(const pen_parserSection& s,
	const T& value,
	const std::string& actualSectionPath = "") const {
    
    if(condType == POSITIVE){
      return value >= static_cast<T>(0.0);
    }

    if(condType == NEGATIVE){
      return value < static_cast<T>(0.0);
    }
    
    if(valType == INT){
      return checkT(value, condInt);
    }else if(valType == DOUBLE){
      return checkT(value, condDouble);
    }else{
      //The value is in a path, read it

      std::string path;
      if(condPath.empty())
	return true;
      if(condPath[0] == '/'){
	//It is an absolute path
	path = condPath;
      }else{
	//It is a relative path
	std::string prefix = actualSectionPath;
	if(!prefix.empty())
	  if(prefix.back() != '/')
	    prefix.append(1,'/');
	path = prefix + condPath;
      }
            
      int auxInt;
      double auxDouble;
      if(s.read(path, auxInt) == INTDATA_SUCCESS){
	return checkT(value, auxInt);
      }else if(s.read(path, auxDouble) == INTDATA_SUCCESS){
	return checkT(value, auxDouble);
      }else{
	return false;
      }
    }
  }

  inline std::string stringify() const {
    std::string result("value must be ");

    switch(condType){
    case POSITIVE:      result.append("positive"); return result;
    case NEGATIVE:      result.append("negative"); return result;
    case LESSER:        result.append("lesser than "); break;
    case LESSER_EQUAL:  result.append("lesser or equal to "); break;
    case GREATER:       result.append("greater than "); break;
    case GREATER_EQUAL: result.append("greater or equal to "); break;
    case NOT_EQUAL:     result.append("not equal to "); break;
    case EQUAL:         result.append("equal to "); break;
    case NO_CONDITION:  return "value has no condition";
    default: return "value has no condition";
    }

    if(valType == INT){
      result.append(std::to_string(condInt));
    }else if(valType == DOUBLE){
      result.append(std::to_string(condDouble));
    }else{
      result.append("data stored in ");
      result.append(condPath);
    }
    
    return result;
  }
  inline std::string stringify(const pen_parserSection& s,
			       const std::string& actualSectionPath = "") const {

    std::string result = stringify();

    std::string path;
    if(condPath.empty())
      return result;
    
    if(condPath[0] == '/'){
      //It is an absolute path
      path = condPath;
    }else{
      //It is a relative path
      std::string prefix = actualSectionPath;
      if(!prefix.empty())
	if(prefix.back() != '/')
	  prefix.append(1,'/');
      path = prefix + condPath;
    }
    
    //Append value from the section
    if(valType == PATH){
      //Read element
      pen_parserElement e;
      if(s.read(path.c_str(),e) != INTDATA_SUCCESS){
	result.append(" (Not defined)");
      }else{
	result.append(" (");
	result.append(e.stringify());
	result.append(")");
      }
    }

    return result;
  }

  inline void clear(){
    condType = NO_CONDITION;
    valType = NONE;
    description.clear();
    condPath.clear();
  }
};

struct pen_readerElement{

  static const std::string valueKey;
  static const std::string descriptionKey;
  static const std::string requiredKey;
  static const std::string requiredTypeKey;
  static const std::string conditionKey;
  
  std::string description;
  pen_parserElement defaultValue;
  pen_readerRequired required;
  std::vector<pen_readerCondition> conditions;

  pen_readerElement() = default;

  //Parse function
  int parse(const pen_parserSection& s,
	    std::string& errorString,
	    const std::string& prefixIn = "");

  inline bool isRequired(const pen_parserSection& s,
			 const std::string& actualSectionPath = "") const{
    return required.required(s,actualSectionPath);
  }

  template<class T>
  inline bool checkConditions(const pen_parserSection& s,
			      const T& value,
			      const std::string& actualSectionPath = "",
			      int* failed = nullptr) const{
    //Iterate over all conditions
    int count = 0;
    for(const pen_readerCondition& c : conditions){
      if(!c.check(s, value, actualSectionPath)){
	*failed = count;
	return false;
      }
      ++count;
    }
    *failed = -1;
    return true;
  }

  std::string stringify(const size_t nSpaces = 0) const;
  inline std::string stringifyExample() const{
    if(defaultValue.readTag() == pen_parserElement::STRING)
      return "\"" + defaultValue.stringify() + "\"";
    else
      return defaultValue.stringify();
  }

  inline void clear(){
    description.clear();
    defaultValue.clear();
    required.clear();
    conditions.clear();
  }

  inline pen_parserElement::types readTag() const{
    return defaultValue.readTag();
  }

  inline const pen_readerCondition& readCondition(const size_t i) const {
    return conditions[i];
  }

  inline const pen_parserElement& readDefault() const {
    return defaultValue;
  }  
};

class pen_readerStorage{

public:

  friend class pen_readerSection;
  
  enum errors{
    SUCCESS = 0,
    UNHANDLED = 1,
  };

  int errorCode;

protected:
  virtual inline void beginRead(){  }
  virtual inline void endRead(){  }
  
  virtual inline int beginSectionFamily(const std::string& /*pathInSection*/,
				 const size_t /*size*/,
				 const unsigned /*verbose*/) { return UNHANDLED; }
  virtual inline int endSectionFamily(const unsigned /*verbose*/) { return UNHANDLED; };
  virtual inline int beginSection(const std::string& /*name*/,
			   const unsigned /*verbose*/) { return UNHANDLED; };
  virtual inline int endSection(const unsigned /*verbose*/) { return UNHANDLED; }
  virtual inline int beginArray(const std::string& /*pathInSection*/,
			 const size_t /*size*/,
			 const unsigned /*verbose*/) { return UNHANDLED; }
  virtual inline int endArray(const unsigned /*verbose*/) { return UNHANDLED; }
  
  virtual inline int storeElement(const std::string& /*pathInSection*/,
			   const pen_parserData& /*element*/,
			   const unsigned /*verbose*/) { return UNHANDLED; }
  virtual inline int storeArrayElement(const std::string& /*pathInSection*/,
				const pen_parserData& /*element*/,
				const size_t /*pos*/,
				const unsigned /*verbose*/) { return UNHANDLED; }
  virtual inline int storeString(const std::string& /*pathInSection*/,
			  const std::string& /*element*/,
			  const unsigned /*verbose*/) { return UNHANDLED; }
  
};

class pen_readerSection{

public:
  static const std::string subsectionArrayKey;
  static const pen_readerElement nullElement;
private:

  bool singleValueSection;
  std::string description;
  std::map<std::string, pen_readerElement> elements;
  std::map<std::string, pen_readerSection> subsections;
  pen_readerRequired required;

  int readSectionNoChecks(const pen_parserSection& in,
			  pen_readerStorage& reader,
			  std::string& errorElement,
			  const std::string& actualSectionPath = "",
			  const unsigned verbose = 2) const;
  
public:

  pen_readerSection() = default;

  int setFormat(const std::string& formatIn,
		std::string& errorString,
		unsigned long& errorLine);
  
  std::string stringify(const size_t nSpaces = 0, const bool printExample = true) const;
  std::string stringifyExample(const std::string& prefixIn) const;
  
  int parse(const pen_parserSection& config,
	    std::string& errorString);

  inline bool exists(const std::string& path) const {

    if(elements.count(path) > 0){
      return true;
    }

    //Check subsections
    for(const auto& s : subsections ){
      std::string subsectionCompletePath = path + subsectionArrayKey;
      if(s.first.find(subsectionCompletePath) == 0){
	//Subsection prefix found, search inside it
	std::string subPath = path.substr(subsectionCompletePath.size());
	return s.second.exists(subPath);
      }
    }

    //Path not found either in elements or subsections
    return false;
  }
  
  inline bool isRequired(const pen_parserSection& s,
			 const std::string& actualSectionPath = "") const{
    if(singleValueSection && elements.size() > 0){
      return elements.cbegin()->second.isRequired(s,actualSectionPath);
    }else{
      return required.required(s,actualSectionPath);
    }
  }
  
  inline const pen_readerElement& readElement(const std::string& path) const{

    auto it = elements.find(path); 
    if(it != elements.end()){
      return it->second;
    }

    //Check subsections
    for(const auto& s : subsections ){
      std::string subsectionCompletePath = path + subsectionArrayKey;
      if(s.first.find(subsectionCompletePath) == 0){
	//Subsection prefix found, search inside it
	std::string subPath = path.substr(subsectionCompletePath.size());
	return s.second.readElement(subPath);
      }
    }

    //Path not found either in elements or subsections
    return nullElement;
  }
    
  inline static constexpr bool isNull(const pen_readerElement& in) {
    return &nullElement == &in;
  }

  bool checkRequired(const pen_parserSection& in,
		     std::string& missingElement,
		     const std::string& actualSectionPath = "") const;  

  bool checkSection(const pen_parserSection& in,
		    std::string& errorElement,
		    std::string& errorElementSpecs,
		    const std::string& actualSectionPath = "") const;
  
  int read(const pen_parserSection& in,
	   pen_readerStorage& reader,
	   std::string& errorElement,
	   std::string& errorElementSpecs,
	   const unsigned verbose) const;
  
  inline void clear(){
    singleValueSection = false;
    description.clear();
    elements.clear();
    subsections.clear();
  }
  
  // ** Static methods

  //Create reader sections form format string
  inline static pen_readerSection
  createFormat(const std::string& formatIn,
	       std::string& errorString,
	       unsigned long& errorLine,
	       int& errCode){

    //Parses and returs a reader section with the provided format
    
    pen_readerSection rs;
    errCode = rs.setFormat(formatIn,
			   errorString,
			   errorLine);
    return rs;
  }
  
  inline static pen_readerSection
  createFormat(const std::string& formatIn){

    //Parses and returs a reader section with the provided format,
    //but throws an exception on fail
    
    std::string errorString;
    unsigned long errorLine;
    int errorCode = 0;
    pen_readerSection rs = createFormat(formatIn,
					errorString, errorLine,
					errorCode);

    if(errorCode != INTDATA_SUCCESS){
      throw std::runtime_error("Error creating format:\n"
			       "Error code: " + std::to_string(errorCode) + "\n"
			       + (errorLine > 0 ? ("Error line: " +
						   std::to_string(errorLine) +
						   "\n") : "") +
			       "Error element: " + errorString + "\n"
			       "Error message: " + pen_parserError(errorCode) + "\n"
			       "Provided format: \n" + formatIn);
    }
    return rs;
  }

  // * Static methods for configurable objects

  //Tests if a object is configurable
  template<class configurableObj>
  static constexpr bool isConfigurable(){
    return pen_format<configurableObj>::format != nullptr;
  }
  
  //Read a pen parser section for a specific object
  template<class configurableObj>
  static const pen_readerSection& readObjectSection(){
    static_assert(isConfigurable<configurableObj>(),
		  "readReaderSection: Error: The provided"
		  " type is not configurable i.e. 'pen_format'"
		  " has not been defined.");

    //Create the reader section with the corresponding format
    static const pen_readerSection rs =
      createFormat(pen_format<configurableObj>::format);

    //Return the reader section
    return rs;
  }
  
  template<class configurableObj>
  static inline int read(const pen_parserSection& in,
			 pen_readerStorage& reader,
			 std::string& errorElement,
			 std::string& errorElementSpecs,
			 const unsigned verbose){

    return readObjectSection<configurableObj>().read(in,
						     reader,
						     errorElement,
						     errorElementSpecs,
						     verbose);
    
  }
  
};

//Define a auxiliary classes with the common functionality of each configurable object

template<class configurableObj>
struct pen_configReader : public pen_readerStorage{
  
  static constexpr bool isConfigurable(){
    return pen_readerSection::isConfigurable<configurableObj>();
  }

  inline int read(const pen_parserSection& config,
		  const unsigned verbose){
    
    std::string errorElement, errorElementSpecs;
    return pen_readerSection::read<configurableObj>(config, *this,
						    errorElement,
						    errorElementSpecs,
						    verbose);
  }
  
};

#endif
