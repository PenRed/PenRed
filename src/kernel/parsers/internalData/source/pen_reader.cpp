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


#include "pen_reader.hh"


// ** Reader classes and functions


const std::string pen_readerRequired::typeKey = "type";
const std::string pen_readerRequired::valueKey = "value";

int pen_readerRequired::parse(const pen_parserSection& s, const std::string& prefixIn){

  std::string prefix = prefixIn;
  if(!prefix.empty()){
    if(prefix.back() != '/'){
      //Append a slash
      prefix.append(1,'/');
    }
  }
  
  //Try to read required type
  std::string reqTypeString;
  if(s.read((prefix + typeKey).c_str(), reqTypeString) != INTDATA_SUCCESS){
    return INTDATA_READER_REQUIRED_WITH_NO_TYPE;
  }

  //Convert to lower
  std::transform(reqTypeString.begin(), reqTypeString.end(), reqTypeString.begin(),
		 [](unsigned char c){ return std::tolower(c); });

  //Check type
  reqType = REQUIRED;
  if(reqTypeString.compare("required") == 0){
    reqType = REQUIRED;
  }else if(reqTypeString.compare("required_if") == 0){
    reqType = REQUIRED_IF;
  }else if(reqTypeString.compare("required_if_exist") == 0){
    reqType = REQUIRED_IF_EXIST;
  }else if(reqTypeString.compare("optional") == 0){
    reqType = OPTIONAL;
  }else if(reqTypeString.compare("optional_if") == 0){
    reqType = OPTIONAL_IF;
  }else if(reqTypeString.compare("optional_if_exist") == 0){
    reqType = OPTIONAL_IF_EXIST;
  }else{
    //Unknown type
    return INTDATA_READER_REQUIRED_UNKNOWN_TYPE;
  }

  //Check if a path is required for this condition
  if(reqType == REQUIRED ||
     reqType == OPTIONAL){
    //No path required
    reqPath.clear();
    return INTDATA_SUCCESS;
  }

  //Read the test path
  if(s.read((prefix + valueKey).c_str(), reqPath) != INTDATA_SUCCESS){
    //Invalid value
    reqPath.clear();
    return INTDATA_READER_REQUIRED_INVALID_VALUE;
  }
  
  return INTDATA_SUCCESS;
}

bool pen_readerRequired::required(const pen_parserSection& s,
				  const std::string& actualSectionPath) const{

  if(reqType == REQUIRED){
    return true;
  }

  if(reqType == OPTIONAL){
    return false;
  }

  //Construct the required path
  std::string finalReqPath;
  if(reqPath[0] != '/'){
    //It is a relative path
    if(!actualSectionPath.empty()){
      finalReqPath = actualSectionPath;
      if(actualSectionPath.back() != '/')
	finalReqPath.append(1,'/');
    }
    finalReqPath += reqPath;
  }else{
    //It is an absolute path, no modifications needed
    finalReqPath = reqPath;
  }

  //Check if the path exists
  bool exists = s.exists(finalReqPath);

  if(reqType == REQUIRED_IF_EXIST)
    return  exists;
  if(reqType == OPTIONAL_IF_EXIST)
    return !exists;
    
  //Try to get a bolean value from stored path
  bool value;
  if(s.read(finalReqPath.c_str(),value) != INTDATA_SUCCESS){
    value = false;
  }
    
  if(reqType == REQUIRED_IF)
    return value;
  if(reqType == OPTIONAL_IF)
    return !value;

  //Required by default
  return true;
}

const std::string pen_readerCondition::typeKey = "type";
const std::string pen_readerCondition::valueKey = "value";
const std::string pen_readerCondition::descriptionKey = "description";

int pen_readerCondition::parse(const pen_parserSection& s,
			       std::string& errorString,
			       const std::string& prefixIn){

  std::string prefix = prefixIn;
  if(!prefix.empty()){
    if(prefix.back() != '/'){
      //Append a slash
      prefix.append(1,'/');
    }
  }
  
  //Try to read condition type
  std::string condTypeString;
  if(s.read((prefix + typeKey).c_str(), condTypeString) != INTDATA_SUCCESS){
    errorString = prefix + typeKey;
    return INTDATA_READER_CONDITION_WITH_NO_TYPE;
  }

  //Convert to lower
  std::transform(condTypeString.begin(), condTypeString.end(), condTypeString.begin(),
		 [](unsigned char c){ return std::tolower(c); });

  //Check type
  condType = NO_CONDITION;
  if(condTypeString.compare("positive") == 0){
    condType = POSITIVE;
  }else if(condTypeString.compare("negative") == 0){
    condType = NEGATIVE;
  }else if(condTypeString.compare("lesser") == 0){
    condType = LESSER;
  }else if(condTypeString.compare("greater") == 0){
    condType = GREATER;
  }else if(condTypeString.compare("not_equal") == 0){
    condType = NOT_EQUAL;
  }else if(condTypeString.compare("equal") == 0){
    condType = EQUAL;
  }else{
    //Unknown type
    errorString = prefix + typeKey + " " + condTypeString;
    return INTDATA_READER_CONDITION_UNKNOWN_TYPE;
  }

  //Check if a description is provided
  if(s.read((prefix + descriptionKey).c_str(), description) != INTDATA_SUCCESS){
    description.clear();
  }

  // ** Non value required conditions
  if(condType == POSITIVE ||
     condType == NEGATIVE){
    //No value required for this condition
    valType = NONE;
    return INTDATA_SUCCESS;    
  }

  // ** Value based conditions (LESSER, GREATER, NOT_EQUAL, EQUAL)
  pen_parserData data;
  if(s.read((prefix + valueKey).c_str(), data) != INTDATA_SUCCESS){
    //Is not a raw value, must be a string where a path is specified
    if(s.read((prefix + valueKey).c_str(), condPath) != INTDATA_SUCCESS){
      //Invalid value
      condPath.clear();
      errorString = prefix + valueKey;
      return INTDATA_READER_CONDITION_INVALID_VALUE;
    }
    valType = PATH;
      
  }else{
    //A constant value is provided
    condPath.clear();
    //Check if it is an integer or a double
    if(data.readTag() == pen_parserData::DOUBLE){
      valType = DOUBLE;
      data.read(condDouble);
    }else if(data.readTag() == pen_parserData::INT){
      valType = INT;
      data.read(condInt);	
    }else{
      errorString = prefix + data.stringify();      
      return INTDATA_READER_CONDITION_INVALID_VALUE;
    }
  }
  
  return INTDATA_SUCCESS;
}


const std::string pen_readerElement::valueKey       = "reader-value";
const std::string pen_readerElement::descriptionKey = "reader-description";
const std::string pen_readerElement::requiredKey    = "reader-required";
const std::string pen_readerElement::requiredTypeKey= "reader-required/type";
const std::string pen_readerElement::conditionKey   = "reader-conditions";

int pen_readerElement::parse(const pen_parserSection& s,
			     std::string& errorString,
			     const std::string& prefixIn){
  
  //Clear
  clear();
  
  std::string prefix = prefixIn;
  if(!prefix.empty()){
    if(prefix.back() != '/'){
      //Append a slash
      prefix.append(1,'/');
    }
  }
    
  // ** Value section
  //Try to read default value
  if(s.read((prefix + valueKey).c_str(), defaultValue) != INTDATA_SUCCESS){
    //Unable to read default value. It is not a element
    errorString = prefix + valueKey;
    return INTDATA_READER_NOT_A_ELEMENT;
  }

  // ** Description section
  //It is a element, try to read the description
  if(s.read((prefix + descriptionKey).c_str(), description) != INTDATA_SUCCESS){
    //Description not provided
    description.assign("");
  }

  // ** "Required" section
  //Try to read and parse the "required" section
  pen_parserSection reqSec;
  if(s.readSubsection((prefix + requiredKey).c_str(),reqSec) == INTDATA_SUCCESS){
    int err = required.parse(reqSec);
    errorString = prefix + requiredKey;
    if(err != INTDATA_SUCCESS)
      return err;
  }

  // ** Conditions section
  //Clear previous conditions
  conditions.clear();
  
  //Read conditions
  std::vector<std::string> conditionKeys;
  s.ls((prefix + conditionKey).c_str(),conditionKeys);
  if(conditionKeys.size() > 0){

    //Resize conditions vector
    conditions.resize(conditionKeys.size());
    
    for(size_t i = 0; i < conditionKeys.size(); ++i){
      
      std::string condKey = prefix + pen_readerElement::conditionKey + "/" + conditionKeys[i];
      
      //Ensure that all condition keys are subsections
      pen_parserSection condSection;
      if(s.read(condKey,condSection) != INTDATA_SUCCESS){
	conditions.clear();
	errorString = condKey;
	return INTDATA_READER_CONDITION_NOT_A_SECTION;
      }

      //Parse section
      int err = conditions[i].parse(condSection, errorString);
      if(err != INTDATA_SUCCESS){
	//Unable to parse condition
	conditions.clear();
	errorString = condKey + ":\n";
	errorString = condSection.stringify();
	return err;
      }
    }
  }
  
  return INTDATA_SUCCESS;
}

template<> inline bool
pen_readerElement::checkConditions(const pen_parserSection& s,
				   const pen_parserData& d,
				   const std::string& actualSectionPath,
				   int* failed) const {
  //Check the type
  switch(d.readTag()){
  case pen_parserData::INT: {
    int aux = d;
    return checkConditions(s,aux,actualSectionPath,failed);
  }
  case pen_parserData::DOUBLE: {
    double aux = d;
    return checkConditions(s,aux,actualSectionPath,failed);
  }
  case pen_parserData::BOOL: {
    bool aux = d;
    return checkConditions(s,aux,actualSectionPath,failed);
  }
  case pen_parserData::CHAR: {
    char aux = d;
    return checkConditions(s,aux,actualSectionPath,failed);
  }
  default: return false;
  }
}
  
template<> inline bool
pen_readerElement::checkConditions(const pen_parserSection& s,
				   const pen_parserElement& e,
				   const std::string& actualSectionPath,
				   int* failed) const {

  //Check the type
  if(e.readTag() != pen_parserElement::SCALAR)
    return true; //Actually, only conditions to scalar types can be applied
    
    
  pen_parserData data;
  if(e.read(data) != INTDATA_SUCCESS)
    return false;
  else
    return checkConditions(s,data,actualSectionPath,failed);
}

std::string pen_readerElement::stringify(const size_t nSpaces) const {

  std::string result;
  std::string beginLine(nSpaces, ' ');
  beginLine = "# " + beginLine;

  //Process description
  if(!description.empty()){
    result += "#\n";

    //Get next endLine in description string
    std::string::size_type endLinePos = description.find("\\n");
    result += beginLine + description.substr(0,endLinePos) + "\n";
    while(endLinePos != std::string::npos){
      std::string::size_type initPos = endLinePos;
      endLinePos = description.find("\\n",initPos+1);
      result += beginLine + description.substr(initPos+2,endLinePos) + "\n";
    }
  }

  //Process required condition
  result += "#\n" + beginLine + required.stringify() + "\n";

  //Process conditions
  if(conditions.size() > 0){
    result += "#\n" + beginLine + "Must satisfy the following conditions:\n#";
    for(const pen_readerCondition& cond : conditions){
      result += "\n" + beginLine + "  - " + cond.stringify();
    }
    result += "\n#\n";
  }else{
    result += "#\n";
  }
    
  //Process value
  result += beginLine + "Default value: " + defaultValue.stringify();

  return result;
}

const std::string pen_readerSection::subsectionArrayKey = "/${subsection}/";

int pen_readerSection::parse(const pen_parserSection& config,
			     std::string& errorString){

  //Clear previous configuration
  clear();
  
  //Iterate over config elements
  const pen_parserSection::elementMap& configElements = config.readMap();
  //Get initial iterator
  auto it = configElements.cbegin();
  //Iterate until map end
  while(it != configElements.cend()){

    //Get the corresponding path
    const std::string& path = it->first;

    //Check if is the section description
    if(path.compare(pen_readerElement::descriptionKey) == 0){
      //Read the description key
      description = it->second.stringify();
      ++it;
      continue;
    }
    
    //Check if is the section 'required'
    if(path.compare(pen_readerElement::requiredTypeKey) == 0){

      //Parse the 'required' section
      required.parse(config, pen_readerElement::requiredKey);
      ++it;
      continue;
    }

    //Check if is the 'value' key
    if(path.compare(pen_readerElement::valueKey) == 0){
      //The section contains no deeper keys, defines directly
      //a value assigned to the section name

      //Ensure this section contains not more elements
      if(elements.size() > 0){
	errorString = path;
	return INTDATA_READER_SINGLE_VALUE_SECTION_WITH_MULTIPLE_KEYS;
      }
      
      singleValueSection = true;
      int err = elements[""].parse(config,
				   errorString,
				   "");
      if(err != INTDATA_SUCCESS){
	return err;
      }
      ++it;
      continue;
    }
    

    //Get the position of the last slash
    std::string::size_type plslash = path.rfind("/");
    if(plslash == std::string::npos){
      errorString = path;
      return INTDATA_READER_NOT_A_SECTION;
    }

    //Check if the last key is a reader element value key
    if(path.compare(plslash+1, std::string::npos,
		    pen_readerElement::valueKey) == 0){

      //This key contains an element value.

      //Check if this section is a single value section
      if(singleValueSection){
	errorString = path;
	return INTDATA_READER_SINGLE_VALUE_SECTION_WITH_MULTIPLE_KEYS;	
      }
      
      //Check if it is inside a subsection array
      std::string::size_type psub = path.find(subsectionArrayKey);
      if(psub == std::string::npos){
	//No subsection contained, increase iterator
	++it;

	//Get the element path
	std::string elementPath = path.substr(0,plslash);
	//Parse element
	int err = elements[elementPath].parse(config,
					      errorString,
					      elementPath);
	if(err != INTDATA_SUCCESS){
	  return err;
	}
	
      }else{
	//This element is a member of a subsection array.

	//Extract the subsection prefix
	const std::string subsectionPath = path.substr(0,psub);
	  
	//Read the subsection
	pen_parserSection subS;
	std::string subCompletePath = subsectionPath + subsectionArrayKey;
	if(config.readSubsection(subCompletePath.c_str(), subS) != INTDATA_SUCCESS){
	  errorString = path;
	  return INTDATA_READER_NOT_A_SECTION;
	}
	  
	//Parse it
	int err = subsections[subsectionPath].parse(subS, errorString);
	if(err != INTDATA_SUCCESS){
	  errorString = subsectionPath + subsectionArrayKey + errorString;
	  return err;
	}

	//Skip the remaining keys belonging this section
	++it;
	while(it != configElements.cend()){
	  //Get the corresponding path
	  const std::string& pathSub = it->first;

	  //Compare with subsection prefix
	  if(pathSub.find(subCompletePath) != 0){
	    //Is not in this section, stop skipping
	    break;
	  }
	  ++it;
	}
      }
    }else{
      //Increase iterator
      ++it;
    }
  }

  return INTDATA_SUCCESS;
}

int pen_readerSection::setFormat(const std::string& formatIn,
				 std::string& errorString,
				 unsigned long& errorLine){

  //Create a stream from the input format string
  std::stringstream iss(formatIn);

  //Create a a section and fill it with the provided format
  pen_parserSection formatSection;
  int err = parseStream(iss, formatSection,
			errorString, errorLine);
  if(err != INTDATA_SUCCESS){
    return err;
  }

  errorLine = 0;
  //Parse format section
  return parse(formatSection, errorString);
}

std::string pen_readerSection::stringify(const size_t nSpaces, const bool printExample) const {

  std::string result;

  std::string beginLine(nSpaces, ' ');
  beginLine = "# " + beginLine;

  if(singleValueSection){
    if(elements.size() > 0){
      result += "#\n#\n" + beginLine + " - It is a single value subsection:\n#\n";
      result += elements.begin()->second.stringify(nSpaces + 6) + "\n#\n";
    }
  }else{  
    if(!description.empty()){
      result += "#\n";

      //Get next endLine in description string
      std::string::size_type endLinePos = description.find("\\n");
      result += beginLine + description.substr(0,endLinePos) + "\n";
      while(endLinePos != std::string::npos){
	std::string::size_type initPos = endLinePos;
	endLinePos = description.find("\\n",initPos+1);
	result += beginLine + description.substr(initPos+2,endLinePos) + "\n";
      }
    }

    //Print the section 'required' configuration
    result += "#\n" + beginLine + required.stringify() + "\n#\n";

    if(elements.size() > 0){
      result += "#\n#\n" + beginLine + " - Contained elements:\n#\n";
      for(const auto& e : elements)
	result += beginLine + "   + " + e.first + ": \n#\n" +
	  e.second.stringify(nSpaces + 6) + "\n#\n";
    }
  }
    
  if(subsections.size() > 0){
    result += "#\n#\n" + beginLine + " - Contained subsections:\n#\n";
    std::string beginLineSub(nSpaces + 6,' ');
    beginLineSub = "# " + beginLineSub;
    for(const auto& s : subsections){
      result += beginLineSub + "* Section family: " + s.first + "\n";
      result += s.second.stringify(nSpaces + 6, false) + "\n";
    }
  }

  if(printExample){
    //Print the configuration example
    result += beginLine + "Configuration example:\n\n";
    result += stringifyExample("");
  }
  
  return result;
}

std::string pen_readerSection::stringifyExample(const std::string& prefixIn) const {
  std::string result;
  std::string actualPath = prefixIn;
  if(!actualPath.empty()){
    if(actualPath.back() != '/')
      actualPath.append(1,'/');
  }

  if(singleValueSection){
    if(elements.size() > 0){
      std::string pathExample = prefixIn;
      if(!pathExample.empty())
	pathExample.pop_back();
      result += pathExample + " " +
	elements.begin()->second.stringifyExample() + "\n";
    }
  }else{
    for(const auto& e : elements){
      result += actualPath + e.first + " " + e.second.stringifyExample() + "\n";
    }
  }

  for(const auto& s : subsections){
    result += s.second.stringifyExample(actualPath + s.first + subsectionArrayKey);
  }
  return result;
}


bool pen_readerSection::checkRequired(const pen_parserSection& in,
				      std::string& missingElement,
				      const std::string& actualSectionPath) const {
  
  //Check if the 'required' conditions are fulfield for
  //all elements in the provided pen_parserSection

  std::string prefix = actualSectionPath;
  if(!prefix.empty())
    if(prefix.back() != '/')
      prefix.append(1,'/');
    

  //Get elements one by one
  for(const auto& e : elements){
    //Check if the element exists
    if(!in.exists((prefix + e.first).c_str())){
      //Is not defined, check if is a required value
      if(e.second.isRequired(in, actualSectionPath)){
	//A required element is not provided, finish
	missingElement = prefix + e.first;
	return false;
      }
    }
  }

  //Check sections
  for(const auto& s : subsections){
    //Check if the subsection array prefix exists
    if(!in.exists((prefix + s.first).c_str())){
      //Is not defined, check if is a required subsection
      if(s.second.isRequired(in, actualSectionPath)){
	//A required subsection array is not provided, finish
	missingElement = prefix + s.first;
	return false;	  
      }
    }else{
      //The subsection array is provided. Check each section name
      std::vector<std::string> subSecNames;
      in.ls(s.first.c_str(), subSecNames);
      for(const std::string& subSecName : subSecNames){
	std::string subSecPath = prefix + s.first + "/" + subSecName;
	bool subReq = s.second.checkRequired(in, missingElement, subSecPath);
	if(!subReq){
	  //An element is missing in subsections
	  return false;
	}
      }
    }
  }
    
  return true;
}

bool pen_readerSection::checkSection(const pen_parserSection& in,
				     std::string& errorElement,
				     std::string& errorElementSpecs,
				     const std::string& actualSectionPath) const {
  
  //Check if the 'required' conditions are fulfield for
  //all elements in the provided pen_parserSection

  std::string prefix = actualSectionPath;
  if(!prefix.empty()){
    if(prefix.back() != '/')
      prefix.append(1,'/');
  }
    

  //Get elements one by one
  for(const auto& e : elements){
    //Check if the element exists
    pen_parserElement eIn;
    if(in.read((prefix + e.first).c_str(), eIn) != INTDATA_SUCCESS){
      //Is not defined, check if is a required value
      if(e.second.isRequired(in, actualSectionPath)){
	//A required element is not provided, finish
	errorElement = prefix + e.first + " " + eIn.stringify();
	errorElementSpecs = e.second.stringify();
	return false;
      }
    } else {
      //The element is defined

      //Check type
      if(eIn.readTag() != e.second.readTag()){
	errorElement = prefix + e.first + " " + eIn.stringify();
	errorElementSpecs = e.second.stringify();	
	return false;
      }
      
      //check conditions
      int failedCond = 0;
      if(!e.second.checkConditions(in,eIn,actualSectionPath,&failedCond)){
	//The conditions have not been fulfilled
	errorElement = prefix + e.first + " " + eIn.stringify();
	errorElementSpecs = e.second.stringify();
	errorElementSpecs += "\n\n# Condition failed: " +
	  e.second.readCondition(failedCond).stringify(in,actualSectionPath);	
	return false;
      }
    }
  }

  //Check sections
  for(const auto& s : subsections){
    //Check if the subsection array prefix exists
    std::string subSectionPath = prefix + s.first;
    if(!in.exists(subSectionPath)){
      //Is not defined, check if is a required subsection
      if(s.second.isRequired(in, actualSectionPath)){
	//A required subsection array is not provided, finish
	errorElement = prefix + s.first;
	errorElementSpecs = s.second.stringify();	
	return false;
      }
    }else{
      //The subsection array is provided. Check each section name
      std::vector<std::string> subSecNames;
      in.ls(subSectionPath.c_str(), subSecNames);
      for(const std::string& subSecName : subSecNames){
	std::string subSecPath = subSectionPath + "/" + subSecName;
	bool subReq = s.second.checkSection(in,
					    errorElement, errorElementSpecs,
					    subSecPath);
	if(!subReq){
	  //An element is missing in subsections
	  return false;
	}
      }
    }
  }
    
  return true;
}

int pen_readerSection::readSectionNoChecks(const pen_parserSection& in,
					   pen_readerStorage& reader,
					   std::string& errorElement,
					   const std::string& actualSectionPath,
					   const unsigned verbose) const {
  
  //Check if the 'required' conditions are fulfield for
  //all elements in the provided pen_parserSection

  std::string prefix = actualSectionPath;
  if(!prefix.empty())
    if(prefix.back() != '/')
      prefix.append(1,'/');

  //Get elements one by one
  for(const auto& e : elements){
    //Check if the element exists
    
    pen_parserElement eIn;
    if(in.read((prefix + e.first).c_str(), eIn) != INTDATA_SUCCESS){
      //Element not provided get default value
      eIn = e.second.readDefault();      
    }
    
    //Check tag and store element
    switch(eIn.readTag()){
    case pen_parserElement::SCALAR:{
      //Read data
      pen_parserData data;
      eIn.read(data);
      //Store data
      int err = reader.storeElement(e.first, data, verbose);
      if(err != 0){
	errorElement = prefix + e.first + " " + eIn.stringify();
	reader.errorCode = err;
	return INTDATA_READER_SPECIFIC_READER_ERROR;
      }
      break;
    }
    case pen_parserElement::STRING:{
      //Read string
      std::string data;
      if(eIn.read(data) != INTDATA_SUCCESS){
	printf("Error: Expected string is not a string. "
	       "Please, report this error.\n");
	errorElement = prefix + e.first + " " + eIn.stringify();
	return INTDATA_NOT_A_STRING;
      }
      //Store it
      int err = reader.storeString(e.first, data, verbose);
      if(err != 0){
	errorElement = prefix + e.first + " " + eIn.stringify();
	reader.errorCode = err;
	return INTDATA_READER_SPECIFIC_READER_ERROR;
      }
      break;
    }
    case pen_parserElement::ARRAY:{
      //Read number of elements in the array
      const size_t arraySize = eIn.size();

      //Flag array start
      int err = reader.beginArray(e.first, arraySize, verbose);
      if(err != 0){
	errorElement = prefix + e.first + " " + eIn.stringify();
	reader.errorCode = err;
	return INTDATA_READER_SPECIFIC_READER_ERROR;
      }
	
      //Save each array position
      for(size_t i = 0; i < arraySize; ++i){
	//Read position 'i'
	pen_parserData data;
	eIn.read(data,i);
	//Store position 'i'
	err = reader.storeArrayElement(e.first, data, i, verbose);
	if(err != 0){
	  errorElement = prefix + e.first + " " + eIn.stringify();
	  reader.errorCode = err;
	  return INTDATA_READER_SPECIFIC_READER_ERROR;
	}
      }

      //Flag array end
      err = reader.endArray(verbose);
      if(err != 0){
	errorElement = prefix + e.first + " " + eIn.stringify();
	reader.errorCode = err;
	return INTDATA_READER_SPECIFIC_READER_ERROR;
      }
	
      break;
    }
    }
    
  }

  //All elements stored, go to subsections
  for(const auto& s : subsections){    
    
    //Check if the subsection array prefix exists
    std::string subSectionPath = prefix + s.first;
    if(in.exists(subSectionPath)){
      //The subsection array is provided.

      //Get each section name
      std::vector<std::string> subSecNames;
      in.ls(subSectionPath.c_str(), subSecNames);      

      //Flag family start
      int err = reader.beginSectionFamily(s.first, subSecNames.size(), verbose);
      if(err != 0){
	errorElement = subSectionPath;
	reader.errorCode = err;
	return INTDATA_READER_SPECIFIC_READER_ERROR;
      }
      
      for(const std::string& subSecName : subSecNames){
	std::string subSecPath = subSectionPath + "/" + subSecName;
	
	//Flag subsection start
	err = reader.beginSection(subSecName, verbose);
	if(err != 0){
	  errorElement = prefix + s.first;
	  reader.errorCode = err;
	  return INTDATA_READER_SPECIFIC_READER_ERROR;
	}
	
	int subRet = s.second.readSectionNoChecks(in, reader, errorElement,
						  subSecPath, verbose);
	if(subRet != INTDATA_SUCCESS){
	  //An element read failed in subsections
	  if(verbose > 0){
	    printf("Error reading section '%s'.\n",subSecPath.c_str());
	  }
	  return subRet;
	}

	//Flag subsection end
	err = reader.endSection(verbose);
	if(err != 0){
	  errorElement = prefix + s.first;
	  reader.errorCode = err;
	  return INTDATA_READER_SPECIFIC_READER_ERROR;
	}
      }

      //Flag family end
      err = reader.endSectionFamily(verbose);
      if(err != 0){
	errorElement = prefix + s.first;
	reader.errorCode = err;
	return INTDATA_READER_SPECIFIC_READER_ERROR;
      }
    }
  }
    
  return INTDATA_SUCCESS;
}


int pen_readerSection::read(const pen_parserSection& in,
			    pen_readerStorage& reader,
			    std::string& errorElement,
			    std::string& errorElementSpecs,
			    const unsigned verbose) const{

  errorElement.clear();
  errorElementSpecs.clear();
  
  //Check requirements and conditions of all elements
  if(!checkSection(in, errorElement, errorElementSpecs)){
    if(verbose > 0){
      printf("Error on elements check.\n"
	     "  - Failed element: '%s'\n"
	     "  - Element format:\n%s\n",
	     errorElement.c_str(),
	     errorElementSpecs.c_str());
    }
    return INTDATA_READER_REQUIREMENTS_AND_CONDITIONS_NOT_FULFILLED;
  }
  
  //All reader requirements are satisfied

  //Flag read begin
  reader.beginRead();

  //read the data and set it
  int err = readSectionNoChecks(in, reader, errorElement, "", verbose);

  if(err != pen_readerStorage::SUCCESS){
    if(verbose > 0){
      printf("Error on read returned by the "
	     "specific class reader.\n"
	     "  - Failed element: '%s'\n"
	     "  - Error code: %d\n",
	     errorElement.c_str(), reader.errorCode);
      return err;
    }
  }

  //Flag read end
  reader.endRead();

  return err;
  
}
