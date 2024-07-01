
//#include "errors.h"

#include <stdio.h>
#include <cstdlib>
#include <string.h>

// C  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
namespace PENERROR_mod		// Initialisation and run-time errors.
{
//       SAVE   ! Saves all items in the module.
  char REASON[128];		// Warning/error message.
  int IRETRN = 0;		// Return code.
  int IERSEC = 0;		// > 0 => Secondary stack overflow.
}
// C  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

void ErrorFunction (int NERROR)
{
  using namespace PENERROR_mod;

  IRETRN = NERROR;
  switch(NERROR){
  case 0: break;
  case 900:
    strcpy (REASON,
	    "Fatal Error reading input: Bad value.");
    break;
  case 901:
    strcpy (REASON,
	    "Fatal Error missing data base file.");
    break;
  case 1001:
    strcpy (REASON,
	    "Fatal Error in PEINIT: Electron absorption energy less than 50 eV.");
    break;
  case 1002:
    strcpy (REASON,
	    "Fatal Error in PEINIT: Photon absorption energy less than 50 eV.");
    break;
  case 1003:
    strcpy (REASON,
	    "Fatal Error in PEINIT: Positron absorption energy less than 50 eV.");
    break;
  case 1004:
    strcpy (REASON,
	    "Fatal Error in PEINIT: The energy interval is too narrow");
    break;
  case 1005:
    strcpy (REASON,
	    "Fatal Error in PEINIT: Too many materials. Edit 'penelope.cpp' and change MAXMAT parameter");
    break;
  case 1006:
    strcpy (REASON,
	    "Fatal Error in PEINIT: The material data file could not be opened");
    break;
  case 1007:
    strcpy (REASON,
	    "Fatal Error in EGRID: The energy interval is too narrow");
    break;
  case 1008:
    strcpy (REASON, "Fatal Error in PEMATR. Too many materials");
    break;
  case 1009:
    strcpy (REASON, "Fatal Error in PEMATR. Corrupt material data file");
    break;
  case 1010:
    strcpy (REASON,
	    "Fatal Error in PEMATR. Too many elements in the material composition. The maximum number of elements is 30");
    break;
  case 1011:
    strcpy (REASON, "Fatal Error in PEMATR. Too many oscillators.");
    break;
  case 1012:
    strcpy (REASON, "Fatal Error in PEMATR. Inconsistent oscillator data.");
    break;
  case 1013:
    strcpy (REASON, "Fatal Error in PEMATR. Too many shells");
    break;
  case 1014:
    strcpy (REASON, "Fatal Error in PEMATR. Too many data points (1).");
    break;
  case 1015:
    strcpy (REASON, "Fatal Error in PEMATR. PHMART. No outer shells?");
    break;
  case 1016:
    strcpy (REASON, "Fatal Error in PEMATR. PHMART. No outer shells?");
    break;
  case 1017:
    strcpy (REASON, "Fatal Error in PEMATR. Too many data points (2).");
    break;
  case 1018:
    strcpy (REASON, "Fatal Error in PEMATR. Corrupt material data file.");
    break;
  case 1101:
    strcpy (REASON,
	    "Fatal Error in PEMATW. This element has been declared twice.");
    break;
  case 1102:
    strcpy (REASON,
	    "Fatal Error in PEMATW. Fractions by weight are too small.");
    break;
  case 1103:
    strcpy (REASON,
	    "Fatal Error in PEMATW. This element has been declared twice.");
    break;
  case 1104:
    strcpy (REASON,
	    "Fatal Error in PEMATW. The allowed material ID numbers are 1-300.");
    break;
  case 1105:
    strcpy (REASON,
	    "Fatal Error in PEMATW. NELEM cannot be larger than 30.");
    break;
  case 1106:
    strcpy (REASON, "Fatal Error in PEMATW. NELEM must be positive.");
    break;
  case 1107:
    strcpy (REASON,
	    "Fatal Error in PEMATW. Abnormal termination of file ''pdcompos.pen''.");
    break;
  case 1108:
    strcpy (REASON, "Fatal Error in PEMATW. Wrong atomic number.");
    break;
  case 1109:
    strcpy (REASON,
	    "Fatal Error in PEMATW. STF (atoms/molecule) must be positive.");
    break;
  case 1110:
    strcpy (REASON,
	    "Fatal Error in PEMATW. Element has been declared twice.");
    break;
  case 1111:
    strcpy (REASON, "Fatal Error in PEMATW. The density must be positive.");
    break;
  case 1112:
    strcpy (REASON, "Fatal Error in PEMATW. Too many shells.");
    break;
  case 1113:
    strcpy (REASON, "Fatal Error in PEMATW. Wrong shell number.");
    break;
  case 1114:
    strcpy (REASON, "Fatal Error in PEMATW. Unbalanced charges (element).");
    break;
  case 1115:
    strcpy (REASON,
	    "Fatal Error in PEMATW. Too many oscillators. The parameter NOM should be increased.");
    break;
  case 1116:
    strcpy (REASON,
	    "Fatal Error in PEMATW. Unbalanced charges (compound).");
    break;
  case 1117:
    strcpy (REASON, "Fatal Error in PEMATW. FP is too large.");
    break;
  case 1118:
    strcpy (REASON,
	    "Fatal Error in PEMATW. Inconsistent oscillator strength data.");
    break;
  case 1119:
    strcpy (REASON,
	    "Fatal Error in PEMATW. Inconsistent oscillator-strength data.");
    break;
  case 1120:
    strcpy (REASON,
	    "Fatal Error in PEMATW. Inconsistent oscillator-strength data (2).");
    break;
  case 1121:
    strcpy (REASON,
	    "Fatal Error in PEMATW. Error in grouping the Compton profiles.");
    break;
  case 1201:
    strcpy (REASON, "Fatal Error in START. E out of range.");
    break;
  case 1202:
    strcpy (REASON, "Fatal Error in KNOCK: Incorrect particle type.");
    break;
  case 1203:
    strcpy (REASON,
	    "Fatal Error in STORES: Not enough storage for secondary positrons.");
    break;
  case 1204:
    strcpy (REASON, "Fatal Error in KNOCKF: Incorrect particle type.");
    break;
  case 1205:
    strcpy (REASON, "Fatal Error in JUMPW: The particle is not a photon.");
    break;
  case 1206:
    strcpy (REASON, "Fatal Error in KNOCKW: The particle is not a photon.");
    break;
  case 1301:
    strcpy (REASON, "Fatal Error in EELa0. Negative total cross section.");
    break;
  case 1302:
    strcpy (REASON, "Fatal Error in EELa0. Inconsistent arguments.");
    break;
  case 1303:
    strcpy (REASON, "Fatal Error in EELaW. Wrong file.");
    break;
  case 1304:
    strcpy (REASON, "Fatal Error in ESIaR. Corrupt material data file.");
    break;
  case 1305:
    strcpy (REASON, "Fatal Error in ESIaR. Too many data points.");
    break;
  case 1306:
    strcpy (REASON, "Fatal Error in ESIaR. Too many shells.");
    break;
  case 1307:
    strcpy (REASON, "Fatal Error in ESIaR. Insufficient memory storage.");
    break;
  case 1308:
    strcpy (REASON, "Fatal Error in ESIaW. Corrupt data file.");
    break;
  case 1309:
    strcpy (REASON, "Fatal Error in ESIaW. Too many shells.");
    break;
  case 1310:
    strcpy (REASON, "Fatal Error in PSIaR. Corrupt material data file.");
    break;
  case 1311:
    strcpy (REASON, "Fatal Error in PSIaR. Too many data points.");
    break;
  case 1312:
    strcpy (REASON, "Fatal Error in PSIaR. Too many shells.");
    break;
  case 1313:
    strcpy (REASON, "Fatal Error in PSIaR. Insufficient memory storage.");
    break;
  case 1314:
    strcpy (REASON, "Fatal Error in PSIaW. Corrupt data file.");
    break;
  case 1315:
    strcpy (REASON, "Fatal Error in PSIaW. Too many shells.");
    break;
  case 1316:
    strcpy (REASON, "Fatal Error in EBRaR. EBRR. Inconsistent format.");
    break;
  case 1317:
    strcpy (REASON, "Fatal Error in EBRW. Corrupt file.");
    break;
  case 1318:
    strcpy (REASON, "Fatal Error in EBRW. Check the bremss database file.");
    break;
  case 1319:
    strcpy (REASON, "Fatal Error in RLMOM. Error code 0.");
    break;
  case 1320:
    strcpy (REASON, "Fatal Error in RLMOM. Error code 1.");
    break;
  case 1321:
    strcpy (REASON, "Fatal Error in RLMOM. Error code 2.");
    break;
  case 1322:
    strcpy (REASON, "Fatal Error in RLMOM. Error code 3.");
    break;
  case 1323:
    strcpy (REASON, "Fatal Error in RLMOM. Error code 4.");
    break;
  case 1324:
    strcpy (REASON, "Fatal Error in BRAR. Inconsistent data.");
    break;
  case 1325:
    strcpy (REASON,
	    "Fatal Error in BRAR. Corrupt data file (pdbrang.p08).");
    break;
  case 1326:
    strcpy (REASON,
	    "Fatal Error in BRAW. Corrupt data file (pdbrang.p08).");
    break;
  case 1327:
    strcpy (REASON, "Fatal Error in GRAaR. RITA initialisation error.");
    break;
  case 1328:
    strcpy (REASON, "Fatal Error in GRAaR. RITA initialisation error.");
    break;
  case 1329:
    strcpy (REASON,
	    "Fatal Error in GRAaR. RITA interpolation error is too large.");
    break;
  case 1330:
    strcpy (REASON, "Fatal Error in GRAaW. Corrupt file.");
    break;
  case 1331:
    strcpy (REASON, "Fatal Error in GRAaW. Corrupt file.");
    break;
  case 1332:
    strcpy (REASON, "Fatal Error in GRAaW. Corrupt file.");
    break;
  case 1333:
    strcpy (REASON, "Fatal Error in GPHaR. Corrupt material data file.");
    break;
  case 1334:
    strcpy (REASON, "Fatal Error in GPHaR. Too many data points.");
    break;
  case 1335:
    strcpy (REASON, "Fatal Error in GPHaR. Too many shells.");
    break;
  case 1336:
    strcpy (REASON, "Fatal Error in GPHaR. Insufficient memory storage.");
    break;
  case 1337:
    strcpy (REASON, "Fatal Error in GPHaW. Too many energies.");
    break;
  case 1338:
    strcpy (REASON, "Fatal Error in GPHaW. Corrupt file.");
    break;
  case 1339:
    strcpy (REASON, "Fatal Error in GPHaW. Too many shells.");
    break;
  case 1340:
    strcpy (REASON, "Fatal Error in GPPaW. Corrupt file.");
    break;
  case 1341:
    strcpy (REASON, "Fatal Error in RELAXR. NTRAN needs to be increased.");
    break;
  case 1342:
    strcpy (REASON, "Fatal Error in RELAXR. Insufficient memory storage.");
    break;
  case 1343:
    strcpy (REASON,
	    "Fatal Error in RELAXR. Negative transition probability?");
    break;
  case 1344:
    strcpy (REASON, "Fatal Error in RELAXR. Rounding error is too large.");
    break;
  case 1345:
    strcpy (REASON, "Fatal Error in RELAXW. The element is not loaded.");
    break;
  case 1346:
    strcpy (REASON, "Fatal Error in RELAXW. NM needs to be increased.");
    break;
  case 1347:
    strcpy (REASON,
	    "Fatal Error in EELdW. Electron cross section data are corrupt.");
    break;
  case 1348:
    strcpy (REASON,
	    "Fatal Error in EELdW. Positron cross section data are corrupt.");
    break;
  case 1349:
    strcpy (REASON,
	    "Fatal Error in EELdR. Error reading electron elastic DCS data.");
    break;
  case 1350:
    strcpy (REASON,
	    "Fatal Error in EELdR. Electron cross section data are corrupt.");
    break;
  case 1351:
    strcpy (REASON,
	    "Fatal Error in EELdR. Error reading positron elastic DCS data.");
    break;
  case 1352:
    strcpy (REASON,
	    "Fatal Error in EELdR. Positron cross section data are corrupt.");
    break;
  case 1353:
    strcpy (REASON, "Fatal Error in EELdR. RITA initialisation error.");
    break;
  case 1354:
    strcpy (REASON, "Fatal Error in EELdR. RITA initialisation error.");
    break;
  case 1355:
    strcpy (REASON, "Fatal Error in ELINIT. Corrupt data file.");
    break;
  case 1356:
    strcpy (REASON, "Fatal Error in ELINIT. Corrupt data file.");
    break;
  case 1357:
    strcpy (REASON, "Fatal Error in DCSEL. Energy out of range.");
    break;
  case 1401:
    strcpy (REASON, "Fatal Error in IRND0. Negative point probability.");
    break;
  case 1402:
    strcpy (REASON, "Fatal Error in RITAI0: N must be larger than 8.");
    break;
  case 1403:
    strcpy (REASON, "Fatal Error in RITAI0: N must be less than NM=512.");
    break;
  case 1404:
    strcpy (REASON,
	    "Fatal Error in RITAI0: XLOW must be larger than XHIGH.");
    break;
  case 1405:
    strcpy (REASON, "Fatal Error in RITAI0: XLOW and NU are negative.");
    break;
  case 1901:
    strcpy (REASON,
	    "Fatal Error in MERGE2. Increase the value of the parameter NP.");
    break;
  case 1902:
    strcpy (REASON,
	    "Fatal Error in MERGE2. Increase the value of the parameter NP.");
    break;
  case 1903:
    strcpy (REASON,
	    "Fatal Error in SORT2. Increase the value of the parameter NP.");
    break;
  case 1904:
    strcpy (REASON, "Fatal Error in SPLINE: N is less than 4.");
    break;
  case 1905:
    strcpy (REASON,
	    "Fatal Error in SPLINE: X values not in increasing order.");
    break;
  case 1906:
    strcpy (REASON, "Fatal Error in SINTEG. Integral limits out of range.");
    break;
  case 1907:
    strcpy (REASON, "Fatal Error int SLAG6: too few data points.");
    break;
  case 1908:
    strcpy (REASON, "Fatal Error in RMOMX. Error code 1.");
    break;
  case 1909:
    strcpy (REASON, "Fatal Error in RMOMX. Error code 2.");
    break;
  case 1910:
    strcpy (REASON, "Fatal Error in RMOMX. Error code 3.");
    break;
  case 1911:
    strcpy (REASON, "Fatal Error in RMOMX. Error code 4.");
    break;
  case 1912:
    strcpy (REASON, "Fatal Error in RNDG30: initialisation error (1).");
    break;
  case 1913:
    strcpy (REASON, "Fatal Error in RNDG30: initialisation error (2).");
    break;
  case 2001:
    strcpy (REASON,
	    "Fatal Error in GEOMIN: The input and output units must be different.");
    break;
  case 2002:
    strcpy (REASON, "Fatal Error in GEOMIN: What do you mean?");
    break;
  case 2003:
    strcpy (REASON, "Fatal Error in GEOMIN: Incorrect label format.");
    break;
  case 2004:
    strcpy (REASON, "Fatal Error in GEOMIN: Same label for two surfaces.");
    break;
  case 2005:
    strcpy (REASON,
	    "Fatal Error in GEOMIN: The parameter NS must be increased.");
    break;
  case 2006:
    strcpy (REASON,
	    "Fatal Error in GEOMIN: Incorrect format of surface indices.");
    break;
  case 2007:
    strcpy (REASON, "Fatal Error in GEOMIN: Incorrect surface indices.");
    break;
  case 2008:
    strcpy (REASON,
	    "Fatal Error in GEOMIN: NPINP is too small (check PARINP).");
    break;
  case 2009:
    strcpy (REASON,
	    "Fatal Error in GEOMIN: Scale factor less than 1.0E-15.)");
    break;
  case 2010:
    strcpy (REASON,
	    "Fatal Error in GEOMIN: Scale factor less than 1.0E-15.)");
    break;
  case 2011:
    strcpy (REASON,
	    "Fatal Error in GEOMIN: Scale factor less than 1.0E-15.)");
    break;
  case 2012:
    strcpy (REASON, "Fatal Error in GEOMIN: What do you mean?");
    break;
  case 2013:
    strcpy (REASON,
	    "Fatal Error in GEOMIN: NPINP is too small (check PARINP)");
    break;
  case 2014:
    strcpy (REASON, "Fatal Error in GEOMIN: What do you mean?");
    break;
  case 2015:
    strcpy (REASON,
	    "Fatal Error in GEOMIN: NPINP is too small (check PARINP)");
    break;
  case 2016:
    strcpy (REASON, "Fatal Error in GEOMIN: What do you mean?");
    break;
  case 2017:
    strcpy (REASON, "Fatal Error in GEOMIN: Incorrect label format.");
    break;
  case 2018:
    strcpy (REASON,
	    "Fatal Error in GEOMIN: Same label for two bodies (or modules).");
    break;
  case 2019:
    strcpy (REASON,
	    "Fatal Error in GEOMIN: The parameter NB must be increased.");
    break;
  case 2020:
    strcpy (REASON,
	    "Fatal Error in GEOMIN: Incorrect material definition line.");
    break;
  case 2021:
    strcpy (REASON, "Fatal Error in GEOMIN: Undefined surface label.");
    break;
  case 2022:
    strcpy (REASON,
	    "Fatal Error in GEOMIN: The last limiting surface has been defined twice.");
    break;
  case 2023:
    strcpy (REASON,
	    "Fatal Error in GEOMIN: The number of limiting surfaces is too large.");
    break;
  case 2024:
    strcpy (REASON, "Fatal Error in GEOMIN: Check side pointer value.");
    break;
  case 2025:
    strcpy (REASON, "Fatal Error in GEOMIN: Undefined body label.");
    break;
  case 2026:
    strcpy (REASON, "Fatal Error in GEOMIN: This body is a module.");
    break;
  case 2027:
    strcpy (REASON,
	    "Fatal Error in GEOMIN: The number of limiting surfaces is too large.");
    break;
  case 2028:
    strcpy (REASON, "Fatal Error in GEOMIN: Undefined body label.");
    break;
  case 2029:
    strcpy (REASON, "Fatal Error in GEOMIN: This module is a body.");
    break;
  case 2030:
    strcpy (REASON,
	    "Fatal Error in GEOMIN: The number of limiting surfaces is too large.");
    break;
  case 2031:
    strcpy (REASON, "Fatal Error in GEOMIN: What do you mean?");
    break;
  case 2032:
    strcpy (REASON, "Fatal Error in GEOMIN: Incorrect label format.");
    break;
  case 2033:
    strcpy (REASON,
	    "Fatal Error in GEOMIN: Same label for two bodies (or modules).");
    break;
  case 2034:
    strcpy (REASON,
	    "Fatal Error in GEOMIN: The parameter NB must be increased.");
    break;
  case 2035:
    strcpy (REASON,
	    "Fatal Error in GEOMIN: Incorrect material definition line.");
    break;
  case 2036:
    strcpy (REASON, "Fatal Error in GEOMIN: Undefined surface label.");
    break;
  case 2037:
    strcpy (REASON,
	    "Fatal Error in GEOMIN: The last limiting surface has been defined twice.");
    break;
  case 2038:
    strcpy (REASON,
	    "Fatal Error in GEOMIN: The number of limiting surfaces is too large.");
    break;
  case 2039:
    strcpy (REASON, "Fatal Error in GEOMIN: Check side pointer value.");
    break;
  case 2040:
    strcpy (REASON, "Fatal Error in GEOMIN: Undefined body label.");
    break;
  case 2041:
    strcpy (REASON, "Fatal Error in GEOMIN: This body is a module.");
    break;
  case 2042:
    strcpy (REASON,
	    "Fatal Error in GEOMIN: You are trying to assign two mothers to the last body.");
    break;
  case 2043:
    strcpy (REASON,
	    "Fatal Error in GEOMIN: The number of limiting surfaces is too large.");
    break;
  case 2044:
    strcpy (REASON, "Fatal Error in GEOMIN: Undefined body label.");
    break;
  case 2045:
    strcpy (REASON, "Fatal Error in GEOMIN: This module is a body.");
    break;
  case 2046:
    strcpy (REASON,
	    "Fatal Error in GEOMIN: You are trying to assign two mothers to the last module.");
    break;
  case 2047:
    strcpy (REASON,
	    "Fatal Error in GEOMIN: The number of limiting surfaces is too large.");
    break;
  case 2048:
    strcpy (REASON, "Fatal Error in GEOMIN: What do you mean?");
    break;
  case 2049:
    strcpy (REASON,
	    "Fatal Error in GEOMIN: NPINP is too small (check PARINP)");
    break;
  case 2050:
    strcpy (REASON, "Fatal Error in GEOMIN: What do you mean?");
    break;
  case 2051:
    strcpy (REASON, "Fatal Error in GEOMIN: Incorrect label format.");
    break;
  case 2052:
    strcpy (REASON,
	    "Fatal Error in GEOMIN: Same label for two bodies or modules.");
    break;
  case 2053:
    strcpy (REASON,
	    "Fatal Error in GEOMIN: The cloned object must be a module.");
    break;
  case 2054:
    strcpy (REASON, "Fatal Error in GEOMIN: This module is not defined.");
    break;
  case 2055:
    strcpy (REASON,
	    "Fatal Error in GEOMIN: The selected object is a body.");
    break;
  case 2056:
    strcpy (REASON,
	    "Fatal Error in GEOMIN: The label does not correspond to a module.");
    break;
  case 2057:
    strcpy (REASON, "Fatal Error in GEOMIN: What do you mean?");
    break;
  case 2058:
    strcpy (REASON,
	    "Fatal Error in GEOMIN: NPINP is too small (check PARINP)");
    break;
  case 2059:
    strcpy (REASON, "Fatal Error in GEOMIN: What do you mean?");
    break;
  case 2060:
    strcpy (REASON,
	    "Fatal Error in GEOMIN: The parameter NS must be increased.");
    break;
  case 2061:
    strcpy (REASON,
	    "Fatal Error in GEOMIN: The parameter NB must be increased.");
    break;
  case 2062:
	
    strcpy (REASON, "Fatal Error in GEOMIN: Something wrong...");
    break;
  case 2063:
	
    strcpy (REASON,
	    "Fatal Error in GEOMIN: The limiting body or module is not yet defined");
    break;
  case 2064:
	
    strcpy (REASON,
	    "Fatal Error in GEOMIN: The limiting body or module is not yet defined");
    break;
  case 2065:
	
    strcpy (REASON, "Fatal Error in GEOMIN: Incorrect label format.");
    break;
  case 2066:
	
    strcpy (REASON, "Fatal Error in GEOMIN: What do you mean?");
    break;
  case 2067:
	
    strcpy (REASON, "Fatal Error in GEOMIN: Too many include levels.");
    break;
  case 2068:
	
    strcpy (REASON, "Fatal Error in GEOMIN: Too many include levels.");
    break;
  case 2069:
	
    strcpy (REASON,
	    "Fatal Error in GEOMIN: The parameter NB must be increased.");
    break;
  case 2070:
	
    strcpy (REASON,
	    "Fatal Error in GEOMIN: The parameter NS must be increased.");
    break;
  case 2071:
	
    strcpy (REASON,
	    "Fatal Error in GEOMIN: The parameter NXG is too small.");
    break;
  case 2072:
	
    strcpy (REASON,
	    "Fatal Error in GEOMIN: Possibly unresolved body or module.");
    break;
  case 2073:
	
    strcpy (REASON, "Fatal Error in GEOMIN: Inconsistent body label.");
    break;
  case 2074:
	
    strcpy (REASON, "Fatal Error in GEOMIN: Inconsistent side pointers.");
    break;
  case 2075:
	
    strcpy (REASON, "Fatal Error in GEOMIN: Wrong input format.");
    break;
    //ERRORS PENEASY
  case 146:
	
    strcpy (REASON,
	    "samplePosition:ERROR: sampling source efficiency is lower than 0.1%%.");
    break;
  case 147:
	
    strcpy (REASON, "iniconfig:ERROR: incorrect section header;");
    break;
  case 148:
	
    strcpy (REASON, "iniconfig:ERROR: too many requested histories");
    break;
  case 149:
	
    strcpy (REASON, "iniconfig:ERROR: refresh interval must be positive.");
    break;
  case 150:
	
    strcpy (REASON, "iniconfig:ERROR: Invalid RNG seeds.");
    break;
  case 151:
	
    strcpy (REASON, "iniconfig:ERROR: unable to open seeds file.");
    break;
  case 152:
	
    strcpy (REASON, "iniconfig:ERROR: unable to open restart file.");
    break;
  case 153:
	
    strcpy (REASON,
	    "iniconfig:ERROR: interval between dumps must be positive.");
    break;
  case 154:
	
    strcpy (REASON, "iniconfig:ERROR: End-Of-Section mark not found");
    break;
  case 155:
	
    strcpy (REASON, "inigeo:ERROR: incorrect section header");
    break;
  case 156:
	
    strcpy (REASON, "inigeo:ERROR: unable to open quadrics file");
    break;
  case 157:
	
    strcpy (REASON,
	    "inigeo:ERROR: too many materials; enlarge MAXMAT parameter and recompile.");
    break;
  case 158:
	
    strcpy (REASON, "inigeo:ERROR: no geometry defined.");
    break;
  case 159:
	
    strcpy (REASON, "inigeo:ERROR: End-Of-Section mark not found");
    break;
  case 160:
	
    strcpy (REASON, "inipen:ERROR: incorrect section header;");
    break;
  case 161:
	
    strcpy (REASON, "inipen:ERROR: Max number of materials exceeded");
    break;
  case 162:
	
    strcpy (REASON, "inipen:ERROR: Invalid MAT index");
    break;
  case 163:
	
    strcpy (REASON, "inipen:ERROR: incomplete list of parameters for MAT");
    break;
  case 164:
	
    strcpy (REASON,
	    "inipen:ERROR: DSMAX must be larger than zero even if electrons are not transported.");
    break;
  case 165:
	
    strcpy (REASON,
	    "inipen:ERROR: There are more materials declared in the geometry file than defined in the config file.");
    break;
  case 166:
	
    strcpy (REASON, "inipen:ERROR: End-Of-Section mark not found");
    break;
  case 167:
	
    strcpy (REASON, "getline:ERROR: unable to read line.");
    break;
  case 168:
	
    strcpy (REASON, "seeki:ERROR: value outside range, xc>x(n):");
    break;
  case 169:
	
    strcpy (REASON, "seeki:ERROR: value outside range, xc<x(1):");
    break;
  case 170:
	
    strcpy (REASON, "BIGSinisrc:ERROR: incorrect section header;");
    break;
  case 171:
	
    strcpy (REASON, "BIGSinisrc:ERROR: Unable to find End-Of-Section");
    break;
  case 172:
	
    strcpy (REASON, "BIGSinisrc:ERROR: expecting to find ON or OFF");
    break;
  case 173:
	
    strcpy (REASON, "BIGSinisrc:ERROR: invalid particle type");
    break;
  case 174:
	
    strcpy (REASON,
	    "BIGSinisrc:ERROR: invalid polarization switch, should be 0 or 1");
    break;
  case 175:
	
    strcpy (REASON,
	    "BIGSinisrc:ERROR: vector P={P1,P2,P3} must be P^2 <= 1");
    break;
  case 176:
	
    strcpy (REASON, "BIGSinisrc:ERROR: null direction.");
    break;
  case 177:
	
    strcpy (REASON, "BIGSinisrc:ERROR: theta1 is less than theta0.");
    break;
  case 178:
	
    strcpy (REASON, "BIGSinisrc:ERROR: invalid interval.");
    break;
  case 179:
	
    strcpy (REASON, "BIGSinisrc:ERROR: Invalid entry. Must be 0 or 1.");
    break;
  case 180:
	
    strcpy (REASON, "BIGSinisrc:ERROR: unable to open spectrum file.");
    break;
  case 181:
	
    strcpy (REASON, "BIGSinisrc:ERROR: invalid entry");
    break;
  case 182:
	
    strcpy (REASON, "BIGSinisrc:ERROR: negative energy");
    break;
  case 183:
	
    strcpy (REASON, "BIGSinisrc:ERROR: decreasing energy");
    break;
  case 184:
	
    strcpy (REASON, "BIGSinisrc:ERROR: too many bins in spectrum;");
    strcpy (REASON, "              enlarge NEMAX");
    break;
  case 185:
	
    strcpy (REASON, "BIGSinisrc:ERROR: at least 1 bin must be defined");
    break;
  case 186:
	
    strcpy (REASON, "BIGSinisrc:ERROR: all probabilities are zero.");
    break;
  case 187:
	
    strcpy (REASON, "BIGSinisrc:ERROR: End-Of-Section mark not found");
    break;
  case 188:
	
    strcpy (REASON,
	    "inisource:ERROR: PSF source ON is incompatible with other sources ON.");
    break;
  case 189:
	
    strcpy (REASON, "EDPinitally:ERROR: incorrect section header;");
    break;
  case 190:
	
    strcpy (REASON,
	    "EDPinitally:ERROR: Unable to find End-Of-Section mark");
    break;
  case 191:
	
    strcpy (REASON, "EDPinitally:ERROR: expecting to find ON or OFF");
    break;
  case 192:
	
    strcpy (REASON, "iniforce:ERROR: incorrect section header");
    break;
  case 193:
	
    strcpy (REASON, "iniforce:ERROR: Unable to find End-Of-Section mark");
    break;
  case 194:
	
    strcpy (REASON, "iniforce:ERROR: expecting to find ON or OFF");
    break;
  case 195:
	
    strcpy (REASON,
	    "iniforce:ERROR: unable to read line containing: MAT,KPAR,ICOL,forcing");
    break;
  case 196:
	
    strcpy (REASON, "iniforce:ERROR: invalid MAT");
    break;
  case 197:
	
    strcpy (REASON, "iniforce:ERROR: proton forcing not implemented yet");
    break;
  case 198:
	
    strcpy (REASON, "iniforce:ERROR: KPAR must be in [1,3]");
    break;
  case 199:
	
    strcpy (REASON, "iniforce:ERROR: ICOL must be in [0,8]");
    break;
  case 200:
	
    strcpy (REASON, "iniforce:ERROR: FORCING must not be < 1");
    break;
  case 201:
	
    strcpy (REASON, "iniforce:ERROR: End-Of-Section mark not found");
    break;
  case 202:
	
    strcpy (REASON, "inisplit:ERROR: incorrect section header");
    break;
  case 203:
	
    strcpy (REASON, "inisplit:ERROR: Unable to find End-Of-Section mark");
    break;
  case 204:
	
    strcpy (REASON, "inisplit:ERROR: expecting to find ON or OFF");
    break;
  case 205:
	
    strcpy (REASON, "inisplit:ERROR: Min weight must be >0.");
    break;
  case 206:
	
    strcpy (REASON,
	    "inisplit:ERROR: Photon polarization is active. Only simple splitting can be used");
    break;
  case 207:
	
    strcpy (REASON, "inisplit:ERROR: Invalid factor");
    break;
  case 208:
	
    strcpy (REASON, "inisplit:ERROR: Invalid sign");
    break;
  case 209:
	
    strcpy (REASON, "inisplit:ERROR: Invalid sector");
    break;
  case 210:
	
    strcpy (REASON, "inisplit:ERROR: End-Of-Section mark not found");
    break;
  case 211:
	
    strcpy (REASON, "splitting:ERROR: internal error");
    break;
  case 212:
	
    strcpy (REASON, "inirussia:ERROR: incorrect section header");
    break;
  case 213:
	
    strcpy (REASON, "inirussia:ERROR: Unable to find End-Of-Section mark");
    break;
  case 214:
	
    strcpy (REASON, "inirussia:ERROR: expecting to find ON or OFF");
    break;
  case 215:
	
    strcpy (REASON, "inirussia:ERROR: Invalid survival probability");
    break;
  case 216:
	
    strcpy (REASON, "inirussia:ERROR: End-Of-Section mark not found");
    break;
  case 217:
	
    strcpy (REASON, "SDDinitially:ERROR: incorrect section header");
    break;
  case 218:
	
    strcpy (REASON, "inivox:ERROR: invalid transparent mat.");
    break;
  case 219:
	
    strcpy (REASON,
	    "inivox:ERROR: there must be ONE body with transparent material.");
    break;
  case 220:
	
    strcpy (REASON,
	    "inivox:ERROR: granularity must be between 2 and maxGranul");
    break;
  case 221:
	
    strcpy (REASON,
	    "inivox:ERROR: internal inconsistency; expecting vacuum");
    break;
  case 222:
	
    strcpy (REASON, "readvox:ERROR: unable to open voxels file");
    break;
  case 223:
	
    strcpy (REASON, "getvoxline:ERROR: unable to read vox file");
    break;
  case 224:
	
    strcpy (REASON, "readvox:ERROR: incorrect section header");
    break;
  case 225:
	
    strcpy (REASON, "readvox:ERROR: invalid no. voxels.");
    break;
  case 226:
	
    strcpy (REASON, "readvox:ERROR: No. voxels exceeds nvoxmax.");
    break;
  case 227:
	
    strcpy (REASON, "readvox:ERROR: not enough memory.");
    break;
  case 228:
	
    strcpy (REASON, "inimassvox:ERROR: not enough memory for temp arrays.");
    break;
  case 229:
	
    strcpy (REASON, "inimassvox:ERROR: not enough memory.");
    break;
  case 230:
	
    strcpy (REASON, "readvox:ERROR: voxel side too small.");
    break;
  case 231:
	
    strcpy (REASON, "readvox:ERROR: VBB too large.");
    break;
  case 232:
	
    strcpy (REASON,
	    "readvox:ERROR: column numbers must be between 1 and maxCol.");
    break;
  case 233:
	
    strcpy (REASON, "readvox:ERROR: End-Of-Section mark not found.");
    break;
  case 234:
	
    strcpy (REASON, "readvox:ERROR: invalid entry at line line.");
    break;
  case 235:
	
    strcpy (REASON, "readvox:ERROR: Line should be blank, line line.");
    break;
  case 236:
	
    strcpy (REASON, "writeMassvox:ERROR: unable to open file to write.");
    break;
  case 237:
	
    strcpy (REASON,
	    "SDDinitally:ERROR: Unable to find End-Of-Section mark.");
    break;
  case 238:
	
    strcpy (REASON, "SDDinitally:ERROR: expecting to find ON or OFF.");
    break;
  case 239:
	
    strcpy (REASON, "SDDinitally:ERROR: invalid entry.");
    break;
  case 240:
	
    strcpy (REASON, "SDDinitally:ERROR: Zero bins defined.");
    break;
  case 241:
	
    strcpy (REASON,
	    "SDDinitally:ERROR: nbinmax max no. of megabins exceeded.");
    break;
  case 242:
	
    strcpy (REASON, "SDDinitally:ERROR: expecting 1,0 or -1.");
    break;
  case 243:
	
    strcpy (REASON, "SDDinitally:ERROR: not enough memory.");
    break;
  case 244:
	
    strcpy (REASON, "SDDinitally:ERROR: End-Of-Section mark not found.");
    break;
  case 245:
	
    strcpy (REASON, "PSFinitally:ERROR: incorrect section header.");
    break;
  case 246:
	
    strcpy (REASON,
	    "PSFinitally:ERROR: Unable to find End-Of-Section mark.");
    break;
  case 247:
	
    strcpy (REASON, "PSFinitally:ERROR: expecting to find ON or OFF.");
    break;
  case 248:
	
    strcpy (REASON,
	    "PSFinitally:ERROR: IAEA PSF format requested but not available.");
    break;
  case 249:
	
    strcpy (REASON, "PSFinitally:ERROR: PSF format must be 0 or 1.");
    break;
  case 250:
	
    strcpy (REASON,
	    "PSFinitally:ERROR: detection material out of range: 1,MAXMAT.");
    break;
  case 251:
	
    strcpy (REASON,
	    "PSFinitally:ERROR: PSF detection material must be a perfect absorbent.");
    break;
  case 252:
	
    strcpy (REASON, "PSFinitally:ERROR: Could not append PSF.");
    break;
  case 253:
	
    strcpy (REASON, "PSFinitally:ERROR: End-Of-Section mark not found.");
    break;
  case 254:
	
    strcpy (REASON, "IAEAiniwrite:ERROR: Unable to open PSF.");
    break;
  case 255:
	
    strcpy (REASON,
	    "IAEAiniwrite:ERROR: Could not open PSF to append; make sure the file exists and it is accessible.");
    break;
  case 256:
	
    strcpy (REASON, "IAEAiniwrite:ERROR: Unable to set an extra variable.");
    break;
  case 257:
	
    strcpy (REASON,
	    "IAEAiniwrite:internalERROR: Extra variable index is out of range.");
    break;
  case 258:
	
    strcpy (REASON,
	    "IAEAiniwrite:internalERROR: Extra variable type is out of range.");
    break;
  case 259:
	
    strcpy (REASON,
	    "IAEAiniwrite:internalERROR: Undefined error while setting an extra variable.");
    break;
  case 260:
	
    strcpy (REASON, "IAEAwrite:ERROR: protons not implemented.");
    break;
  case 261:
	
    strcpy (REASON, "IAEAwrite:ERROR: Invalid KPAR.");
    break;
  case 262:
	
    strcpy (REASON, "IAEAwrite:ERROR: Unable to write particle.");
    break;
  case 263:
	
    strcpy (REASON, "VDDinitially:ERROR: incorrect section header");
    break;
  case 264:
	
    strcpy (REASON,
	    "VDDinitially:ERROR: voxel dose tally is ON but no voxelized geometry has been defined.");
    break;
  case 265:
	
    strcpy (REASON,
	    "VDDinitally:ERROR: Unable to find End-Of-Section mark.");
    break;
  case 266:
	
    strcpy (REASON, "VDDinitally:ERROR: expecting to find ON or OFF.");
    break;
  case 267:
	
    strcpy (REASON, "VDDinitally:ERROR: invalid ROI.");
    break;
  case 268:
	
    strcpy (REASON, "VDDinitally:ERROR: expecting 1,0 or -1.");
    break;
  case 269:
	
    strcpy (REASON, "VDDinitally:ERROR: End-Of-Section mark not found.");
    break;
  case 270:
	
    strcpy (REASON, "VDDinitally:ERROR: not enough memory.");
    break;
  case 271:
	
    strcpy (REASON, "CDDinitally:ERROR: incorrect section header.");
    break;
  case 272:
	
    strcpy (REASON,
	    "CDDinitally:ERROR: Unable to find End-Of-Section mark.");
    break;
  case 273:
	
    strcpy (REASON, "CDDinitally:ERROR: expecting to find ON or OFF.");
    break;
  case 274:
	
    strcpy (REASON, "CDDinitally:ERROR: invalid entry.");
    break;
  case 275:
	
    strcpy (REASON, "CDDinitally:ERROR: Too many bins.");
    break;
  case 276:
	
    strcpy (REASON, "CDDinitally:ERROR: End-Of-Section mark not found.");
    break;
  case 277:
	
    strcpy (REASON, "PSFinisrc:ERROR: incorrect section header.");
    break;
  case 278:
	
    strcpy (REASON, "PSFinisrc:ERROR: Unable to find End-Of-Section mark.");
    break;
  case 279:
	
    strcpy (REASON, "PSFinisrc:ERROR: expecting to find ON or OFF.");
    break;
  case 280:
	
    strcpy (REASON,
	    "PSFinisrc:ERROR: IAEA PSF format requested but not available.");
    break;
  case 281:
	
    strcpy (REASON, "PSFinisrc:ERROR: PSF format must be 0 or 1.");
    break;
  case 282:
	
    strcpy (REASON, "PSFinisrc:ERROR: split < 1.");
    break;
  case 283:
	
    strcpy (REASON,
	    "PSFinisrc:ERROR: invalid entry. VALIDATE field must be 0 or 1.");
    break;
  case 284:
	
    strcpy (REASON, "PSFinisrc:ERROR: cannot open the PSF.");
    break;
  case 285:
	
    strcpy (REASON,
	    "PSFinisrc:ERROR: No. of histories in PSF exceeds nmax.");
    break;
  case 286:
	
    strcpy (REASON, "PSFinisrc:ERROR: invalid KPAR: at line:.");
    break;
  case 287:
	
    strcpy (REASON, "PSFinisrc:ERROR: invalid energy(eV): at line:.");
    break;
  case 288:
	
    strcpy (REASON,
	    "PSFinisrc:ERROR: null vector direction found at line:.");
    break;
  case 289:
	
    strcpy (REASON,
	    "PSFinisrc:ERROR: PSF end-of-file reached, not enough particles to initialize.");
    break;
  case 290:
	
    strcpy (REASON,
	    "PSFinisrc:ERROR: Inconsistency found in PSF incremental history no.");
    break;
  case 291:
	
    strcpy (REASON, "PSFinisrc:ERROR: End-Of-Section mark not found.");
    break;
  case 292:
	
    strcpy (REASON, "getpar:ERROR: unable to read PSF line no.:.");
    break;
  case 293:
	
    strcpy (REASON, "getpar:ERROR: invalid or missing datum in PSF line:.");
    break;
  case 294:
	
    strcpy (REASON, "checkFormat:ERROR: cannot open the PSF.\n");
    break;
  case 295:
	
    strcpy (REASON, "ckeckFormat:ERROR: unable to read first PSF line.\n");
    break;
  case 296:
	
    strcpy (REASON, "checkFormat:ERROR: unable to identify PSF format.\n");
    break;
  case 297:
	
    strcpy (REASON, "getparIAEA:ERROR: unable to read particle no.:.\n");
    break;
  case 298:
	
    strcpy (REASON,
	    "getparIAEA:ERROR: EOF reached when attempting to read particle no.:.\n");
    break;
  case 299:
	
    strcpy (REASON,
	    "getparIAEA:ERROR: undefined error while attempting to read particle no.:.\n");
    break;
  case 300:
	
    strcpy (REASON,
	    "getparIAEA:ERROR: Invalid KPAR: found at particle no.:\n");
    break;
  case 301:
	
    strcpy (REASON, "IAEAiniread:ERROR: unable to open PSF.\n");
    break;
  case 302:
	
    strcpy (REASON,
	    "IAEAiniread:ERROR: Unable to get the total number of particles.\n");
    break;
  case 303:
	
    strcpy (REASON,
	    "IAEAiniread:ERROR: Unable to get the total number of histories.\n");
    break;
  case 304:
	
    strcpy (REASON,
	    "IAEAiniread:ERROR: No. of histories in PSF exceeds nmax.\n");
    break;
  case 305:
	
    strcpy (REASON,
	    "IAEAiniread:ERROR: Unable to get constant variables.\n");
    break;
  case 306:
	
    strcpy (REASON,
	    "IAEAiniread:ERROR: Variable constant index out of range.\n");
    break;
  case 307:
	
    strcpy (REASON,
	    "IAEAiniread:ERROR: More than xtraMax extra int or real variables in PSF.\n");
    break;
  case 308:
	
    strcpy (REASON,
	    "IAEAiniread:ERROR: Unable to check file size and byte order of the PSF.\n");
    break;
  case 309:
	
    strcpy (REASON, "IAEAiniread:ERROR: The function fseek fails.\n");
    break;
  case 310:
	
    strcpy (REASON,
	    "IAEAiniread:ERROR: File size inconsistent with header checksum.\n");
    break;
  case 311:
	
    strcpy (REASON, "IAEAiniread:ERROR: There is a byte order mismatch.\n");
    break;
  case 312:
	
    strcpy (REASON,
	    "IAEAiniread:ERROR: There is a file size and byte order mismatch.\n");
    break;
  case 313:
	
    strcpy (REASON, "IAEAiniread:ERROR: Unidentified error code.\n");
    break;
  case 314:
	
    strcpy (REASON, "IAEAiniread:internalERROR: invalid KPAR.\n");
    break;
  case 315:
	
    strcpy (REASON,
	    "IAEAiniread:ERROR: energy out of range found in particle no.:.\n");
    break;
  case 316:
	
    strcpy (REASON,
	    "IAEAiniread:ERROR: null direction vector found in particle no.:.\n");
    break;
  case 317:
	
    strcpy (REASON,
	    "IAEAiniread:ERROR: No. of particles inferred from PSF differs from that stated in header file.\n");
    break;
  case 318:
	
    strcpy (REASON,
	    "IAEAiniread:ERROR: No. of histories inferred from PSF is greater than that stated in header file.\n");
    break;
  case 319:
	
    strcpy (REASON, "IAEAiniread:ERROR: Unable to get maximum energy.\n");
    break;
  case 320:
	
    strcpy (REASON, "IAEAiniread:ERROR: unable to close PSF.\n");
    break;
  case 321:
	
    strcpy (REASON, "IAEAiniread:ERROR: unable to open PSF.\n");
    break;
  case 322:
	
    strcpy (REASON,
	    "IAEAiniread:ERROR: PSF end-of-file reached, not enough particles to initialize.\n");
    break;
  case 323:
	
    strcpy (REASON,
	    "IAEAiniread:ERROR: Inconsistency found in PSF incremental history.\n");
    break;
  case 324:
	
    strcpy (REASON, "PSFinitally:ERROR: incorrect section header.\n");
    break;
  case 325:
	
    strcpy (REASON,
	    "PSFinitally:ERROR: Unable to find End-Of-Section mark.\n");
    break;
  case 326:
	
    strcpy (REASON, "PSFinitally:ERROR: expecting to find ON or OFF.\n");
    break;
  case 327:
	
    strcpy (REASON,
	    "PSFinitally:ERROR: IAEA PSF format requested but not available.\n");
    break;
  case 328:
	
    strcpy (REASON, "PSFinitally:ERROR: PSF format must be 0 or 1.\n");
    break;
  case 329:
	
    strcpy (REASON,
	    "PSFinitally:ERROR: detection material out of range.\n");
    break;
  case 330:
	
    strcpy (REASON,
	    "PSFinitally:ERROR: PSF detection material must be a perfect absorbent.\n");
    break;
  case 331:
	
    strcpy (REASON, "PSFinitally:ERROR: Could not append PSF.\n");
    break;
  case 332:
	
    strcpy (REASON, "PSFinitally:ERROR: End-Of-Section mark not found.\n");
    break;
  case 333:
	
    strcpy (REASON,
	    "IAEAwriteReport:ERROR: unable to write no. of histories to header file.\n");
    break;
  case 334:
	
    strcpy (REASON,
	    "IAEAwriteReport:ERROR: unable to update header file.\n");
    break;
  case 335:
	
    strcpy (REASON, "FTLinitally:ERROR: incorrect section header.\n");
    break;
  case 336:
	
    strcpy (REASON,
	    "FTLinitally:ERROR: Unable to find End-Of-Section mark.\n");
    break;
  case 337:
	
    strcpy (REASON, "FTLinitally:ERROR: expecting to find ON or OFF.\n");
    break;
  case 338:
	
    strcpy (REASON, "FTLinitally:ERROR: Too many bins.\n");
    break;
  case 339:
	
    strcpy (REASON, "FTLinitally:ERROR: End-Of-Section mark not found.\n");
    break;
  case 340:
	
    strcpy (REASON, "PCSinitally:ERROR: incorrect section header.\n");
    break;
  case 341:
	
    strcpy (REASON,
	    "PCSinitally:ERROR: Unable to find End-Of-Section mark.\n");
    break;
  case 342:
	
    strcpy (REASON, "PCSinitally:ERROR: expecting to find ON or OFF.\n");
    break;
  case 343:
	
    strcpy (REASON, "PCSinitally:ERROR: Too many bins.\n");
    break;
  case 344:
	
    strcpy (REASON, "PCSinitally:ERROR: End-Of-Section mark not found.\n");
    break;
  case 345:
	
    strcpy (REASON, "PTSinitally:ERROR: incorrect section header.\n");
    break;
  case 346:
	
    strcpy (REASON,
	    "PTSinitally:ERROR: Unable to find End-Of-Section mark.\n");
    break;
  case 347:
	
    strcpy (REASON, "PTSinitally:ERROR: expecting to find ON or OFF.\n");
    break;
  case 348:
	
    strcpy (REASON, "PTSinitally:ERROR: cannotopen track data file.\n");
    break;
  case 349:
	
    strcpy (REASON, "PTSinitally:ERROR: End-Of-Section mark not found.\n");
    break;
  case 350:
	
    strcpy (REASON, "SPDinitally:ERROR: incorrect section header.\n");
    break;
  case 351:
	
    strcpy (REASON,
	    "SPDinitally:ERROR: Unable to find End-Of-Section mark.\n");
    break;
  case 352:
	
    strcpy (REASON, "SPDinitally:ERROR: expecting to find ON or OFF.\n");
    break;
  case 353:
	
    strcpy (REASON, "SPDinitally:ERROR: Invalid entry.\n");
    break;
  case 354:
	
    strcpy (REASON, "SPDinitally:ERROR: Too many bins.\n");
    break;
  case 355:
	
    strcpy (REASON, "SPDinitally:ERROR: End-Of-Section mark not found.\n");
    break;
  case 356:
	
    strcpy (REASON, "PHSinitally:ERROR: incorrect section header.\n");
    break;
  case 357:
	
    strcpy (REASON,
	    "PHSinitally:ERROR: Unable to find End-Of-Section mark.\n");
    break;
  case 358:
	
    strcpy (REASON, "PHSinitally:ERROR: expecting to find ON or OFF.\n");
    break;
  case 359:
	
    strcpy (REASON, "PHSinitally:ERROR: Invalid entry.\n");
    break;
  case 360:
	
    strcpy (REASON, "PHSinitally:ERROR: Too many bins.\n");
    break;
  case 361:
	
    strcpy (REASON, "PHSinitally:ERROR: End-Of-Section mark not found.\n");
    break;

  case 362:
	
    strcpy (REASON, "PIDinitally:ERROR: incorrect section header.\n");
    break;
  case 363:
	
    strcpy (REASON,
	    "PIDinitally:ERROR: Unable to find End-Of-Section mark.\n");
    break;
  case 364:
	
    strcpy (REASON, "PIDinitally:ERROR: expecting to find ON or OFF.\n");
    break;
  case 365:
	
    strcpy (REASON, "PIDinitally:ERROR: Invalid entry.\n");
    break;
  case 366:
	
    strcpy (REASON,
	    "PIDinitally:ERROR: detection material must be a perfect absorbent; increase absorption energies above einf.\n");
    break;
  case 367:
	
    strcpy (REASON, "PIDinitally:ERROR: invalid value.\n");
    break;
  case 368:
	
    strcpy (REASON,
	    "PIDinitally:ERROR: pixel size and no. pixels are both zero.\n");
    break;
  case 369:
	
    strcpy (REASON, "PIDinitally:ERROR: max no. megapixels exceeded.\n");
    break;
  case 370:
	
    strcpy (REASON,
	    "PIDinitally:ERROR: emin, emax, nebin invalid values.\n");
    break;
  case 371:
	
    strcpy (REASON, "PIDinitally:ERROR: invalid mode.\n");
    break;
  case 372:
	
    strcpy (REASON,
	    "PIDinitally:ERROR: max no. (Mpixels x E bins) exceeded.\n");
    break;
  case 373:
	
    strcpy (REASON,
	    "PIDinitally:ERROR: c0, c1 parameters cannot be negative.\n");
    break;
  case 374:
	
    strcpy (REASON,
	    "PIDinitally:ERROR: Matrix format is incompatible with energy discriminating mode.\n");
    break;
  case 375:
	
    strcpy (REASON, "PIDinitally:ERROR: invalid format.\n");
    break;
  case 376:
	
    strcpy (REASON, "PIDinitally:ERROR: End-Of-Section mark not found.\n");
    break;
  case 377:
	
    strcpy (REASON, "PIDinitally:ERROR: not enough memory.\n");
    break;
  case 378:
	
    strcpy (REASON,
	    "setDetectFrame:ERROR: detection material not found.\n");
    break;

  case 400:
	
    strcpy (REASON,
	    "PMRDR:ERROR: The input file must begin with the TITLE line.\n");
    break;
  case 401:
	
    strcpy (REASON, "PMRDR:ERROR: Incorrect particle type.\n");
    break;
  case 402:
	
    strcpy (REASON,
	    "PMRDR:ERROR: Source energy spectrum. The number of energy bins is too large.\n");
    break;
  case 403:
	
    strcpy (REASON,
	    "PMRDR:ERROR: The source energy spectrum is not defined.\n");
    break;
  case 404:
	
    strcpy (REASON, "PMRDR:ERROR: The initial energy E0 is too small.\n");
    break;
  case 405:
	
    strcpy (REASON, "PMRDR:ERROR: Incorrect body label.\n");
    break;
  case 406:
	
    strcpy (REASON, "PMRDR:ERROR: THETA must be between 0 and 180 deg.\n");
    break;
  case 407:
	
    strcpy (REASON, "PMRDR:ERROR: PHI must be between 0 and 360 deg.\n");
    break;
  case 408:
	
    strcpy (REASON, "PMRDR:ERROR: ALPHA must be between 0 and 180 deg.\n");
    break;
  case 409:
	
    strcpy (REASON,
	    "PMRDR:ERROR: Inconsistent definition of the primary source.\n");
    break;
  case 410:
	
    strcpy (REASON, "PMRDR:ERROR: Too many phase-space files.\n");
    break;
  case 411:
	
    strcpy (REASON, "PMRDR:ERROR: Inconsistent window end points.\n");
    break;
  case 412:
	
    strcpy (REASON,
	    "PMRDR:ERROR: You have to specify a material file (line MFNAME).\n");
    break;
  case 413:
	
    strcpy (REASON, "PMRDR:ERROR: Wrong number of materials.\n");
    break;
  case 414:
	
    strcpy (REASON, "PMRDR:ERROR: Geometry file could not be opened.\n");
    break;
  case 415:
	
    strcpy (REASON, "PMRDR:ERROR: The PARINP index must be positive.\n");
    break;
  case 416:
	
    strcpy (REASON, "PMRDR:ERROR: Too many modified parameters.\n");
    break;
  case 417:
	
    strcpy (REASON, "PMRDR:ERROR: NMATG must be greater than 0.\n");
    break;
  case 418:
	
    strcpy (REASON, "PMRDR:ERROR: Too many bodies.\n");
    break;
  case 419:
	
    strcpy (REASON, "PMRDR:ERROR: Too many different materials.\n");
    break;
  case 420:
	
    strcpy (REASON, "PMRDR:ERROR: Some source bodies are undefined.\n");
    break;
  case 421:
	
    strcpy (REASON, "PMRDR:ERROR: You have to specify a geometry file.\n");
    break;
  case 422:
	
    strcpy (REASON, "PMRDR:ERROR: Incorrect body number.\n");
    break;
  case 423:
	
    strcpy (REASON, "PMRDR:ERROR: Incorrect KB value.\n");
    break;
  case 424:
	
    strcpy (REASON, "PMRDR:ERROR: Incorrect value of KPAR.\n");
    break;
  case 425:
	
    strcpy (REASON, "PMRDR:ERROR: Incorrect value of ICOL.\n");
    break;
  case 426:
	
    strcpy (REASON, "PMRDR:ERROR: Incorrect weight window limits.\n");
    break;
  case 427:
	
    strcpy (REASON, "PMRDR:ERROR: Incorrect value of IBRSPL.\n");
    break;
  case 428:
	
    strcpy (REASON,
	    "PMRDR:ERROR: Interaction forcing unactive in this body.\n");
    break;
  case 429:
	
    strcpy (REASON, "PMRDR:ERROR: Incorrect value of IXRSPL.\n");
    break;
  case 430:
	
    strcpy (REASON, "PMRDR:ERROR: NBE equal to 0.\n");
    break;
  case 431:
	
    strcpy (REASON, "PMRDR:ERROR: NBTH equal to 0.\n");
    break;
  case 432:
	
    strcpy (REASON, "PMRDR:ERROR: Wrong number of PHI bins.\n");
    break;
  case 433:
	
    strcpy (REASON, "PMRDR:ERROR: Incorrect number of energy bins.\n");
    break;
  case 434:
	
    strcpy (REASON, "PMRDR:ERROR: Incorrect energy limits.\n");
    break;
  case 435:
	
    strcpy (REASON, "PMRDR:ERROR: Wrong IPSF value.\n");
    break;
  case 436:
	
    strcpy (REASON, "PMRDR:ERROR: Wrong IDCUT value.\n");
    break;
  case 437:
	
    strcpy (REASON,
	    "PMRDR:ERROR: No impact detector has been defined yet.\n");
    break;
  case 438:
	
    strcpy (REASON,
	    "PMRDR:ERROR: Only one PSF can be generated in a each run.\n");
    break;
  case 439:
	
    strcpy (REASON, "PMRDR:ERROR: Incorrect number of age bins.\n");
    break;
  case 440:
	
    strcpy (REASON, "PMRDR:ERROR: Incorrect age limits.\n");
    break;
  case 441:
	
    strcpy (REASON, "PMRDR:ERROR: Undefined age distribution limits.\n");
    break;
  case 442:
	
    strcpy (REASON, "PMRDR:ERROR: Incorrect body label.\n");
    break;
  case 443:
	
    strcpy (REASON,
	    "PMRDR:ERROR: A body cannot be part of two detectors.\n");
    break;
  case 444:
	
    strcpy (REASON,
	    "PMRDR:ERROR: A void body cannot be part of a detectors.\n");
    break;
  case 445:
	
    strcpy (REASON, "PMRDR:ERROR: This detector has no active bodies.\n");
    break;
  case 446:
	
    strcpy (REASON, "PMRDR:ERROR: XU must be greater than XL+1.0E-6.\n");
    break;
  case 447:
	
    strcpy (REASON, "PMRDR:ERROR: Incorrect keyword.\n");
    break;
  case 448:
	
    strcpy (REASON, "PMRDR:ERROR: YU must be greater than YL+1.0E-6.\n");
    break;
  case 449:
	
    strcpy (REASON, "PMRDR:ERROR: ZU must be greater than ZL+1.0E-6.\n");
    break;
  case 450:
	
    strcpy (REASON, "PMRDR:ERROR: RU must be greater than 1.0E-6.\n");
    break;
  case 451:
	
    strcpy (REASON,
	    "PMRDR:ERROR: The dump file is corrupted (the TITLE does not match).\n");
    break;
  case 452:
	
    strcpy (REASON, "PMRDR:ERROR: File could not be opened.\n");
    break;
  case 453:
	
    strcpy (REASON, "PMRDR:ERROR: The file is empty or corrupted.\n");
    break;
  case 454:
	
    strcpy (REASON, "ENANG0:ERROR: NBE is too large.\n");
    break;
  case 455:
	
    strcpy (REASON, "ENANG0:ERROR: NBTH is too large.\n");
    break;
  case 456:
	
    strcpy (REASON, "ENANG0:ERROR: NBPH is too large.\n");
    break;
  case 457:
	
    strcpy (REASON,
	    "IMDET0:ERROR: SIMDET: Detector already defined. Detector cannot be defined.\n");
    break;
  case 458:
	
    strcpy (REASON, "IMDET0:ERROR: SIMDET: Too many detectors.\n");
    break;
  case 459:
	
    strcpy (REASON, "IMDET0:ERROR: SIMDET: NB is too large.\n");
    break;
  case 460:
	
    strcpy (REASON,
	    "ENDET0:ERROR: SENDET: Detector already defined. Detector cannot be defined.\n");
    break;
  case 461:
	
    strcpy (REASON, "ENDET0:ERROR: SENDET: Too many detectors.\n");
    break;
  case 462:
	
    strcpy (REASON, "ENDET0:ERROR: SENDET: NB is too large.\n");
    break;
  case 463:
	
    strcpy (REASON, "DOSE0:ERROR: IDOSE should be 1, 2, or 3.\n");
    break;
  case 464:
	
    strcpy (REASON,
	    "DOSE0:ERROR: SDOSE: NBX must be .GT.0. and .LE.NDXM\n");
    break;
  case 465:
	
    strcpy (REASON,
	    "DOSE0:ERROR: SDOSE: NBY must be .GT.0. and .LE.NDYM\n");
    break;
  case 466:
	
    strcpy (REASON,
	    "DOSE0:ERROR: SDOSE: NBZ must be .GT.0. and .LE.NDZM\n");
    break;
  case 1000:
	
    strcpy (REASON, "");
    break;
  default:
    break;
  }
}
