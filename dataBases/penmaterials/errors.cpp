
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

void
ErrorFunction (int NERROR)
{
  using namespace PENERROR_mod;

  IRETRN = NERROR;
  if (NERROR == 0)
    {
      //Acaba correctament
    }
  else if (NERROR == 1001)
    {
      strcpy (REASON,
	      "Fatal Error in PEINIT: Electron absorption energy less than 50 eV.");
    }
  else if (NERROR == 1002)
    {
      strcpy (REASON,
	      "Fatal Error in PEINIT: Photon absorption energy less than 50 eV.");
    }
  else if (NERROR == 1003)
    {
      strcpy (REASON,
	      "Fatal Error in PEINIT: Positron absorption energy less than 50 eV.");
    }
  else if (NERROR == 1004)
    {
      strcpy (REASON,
	      "Fatal Error in PEINIT: The energy interval is too narrow");
    }
  else if (NERROR == 1005)
    {
      strcpy (REASON,
	      "Fatal Error in PEINIT: Too many materials. Edit 'penelope.cpp' and change MAXMAT parameter");
    }
  else if (NERROR == 1006)
    {
      strcpy (REASON,
	      "Fatal Error in PEINIT: The material data file could not be opened");
    }
  else if (NERROR == 1007)
    {
      strcpy (REASON,
	      "Fatal Error in EGRID: The energy interval is too narrow");
    }
  else if (NERROR == 1008)
    {
      strcpy (REASON, "Fatal Error in PEMATR. Too many materials");
    }
  else if (NERROR == 1009)
    {
      strcpy (REASON, "Fatal Error in PEMATR. Corrupt material data file");
    }
  else if (NERROR == 1010)
    {
      strcpy (REASON,
	      "Fatal Error in PEMATR. Too many elements in the material composition. The maximum number of elements is 30");
    }
  else if (NERROR == 1011)
    {
      strcpy (REASON, "Fatal Error in PEMATR. Too many oscillators.");
    }
  else if (NERROR == 1012)
    {
      strcpy (REASON, "Fatal Error in PEMATR. Inconsistent oscillator data.");
    }
  else if (NERROR == 1013)
    {
      strcpy (REASON, "Fatal Error in PEMATR. Too many shells");
    }
  else if (NERROR == 1014)
    {
      strcpy (REASON, "Fatal Error in PEMATR. Too many data points (1).");
    }
  else if (NERROR == 1015)
    {
      strcpy (REASON, "Fatal Error in PEMATR. PHMART. No outer shells?");
    }
  else if (NERROR == 1016)
    {
      strcpy (REASON, "Fatal Error in PEMATR. PHMART. No outer shells?");
    }
  else if (NERROR == 1017)
    {
      strcpy (REASON, "Fatal Error in PEMATR. Too many data points (2).");
    }
  else if (NERROR == 1018)
    {
      strcpy (REASON, "Fatal Error in PEMATR. Corrupt material data file.");
    }
  else if (NERROR == 1101)
    {
      strcpy (REASON,
	      "Fatal Error in PEMATW. This element has been declared twice.");
    }
  else if (NERROR == 1102)
    {
      strcpy (REASON,
	      "Fatal Error in PEMATW. Fractions by weight are too small.");
    }
  else if (NERROR == 1103)
    {
      strcpy (REASON,
	      "Fatal Error in PEMATW. This element has been declared twice.");
    }
  else if (NERROR == 1104)
    {
      strcpy (REASON,
	      "Fatal Error in PEMATW. The allowed material ID numbers are 1-300.");
    }
  else if (NERROR == 1105)
    {
      strcpy (REASON,
	      "Fatal Error in PEMATW. NELEM cannot be larger than 30.");
    }
  else if (NERROR == 1106)
    {
      strcpy (REASON, "Fatal Error in PEMATW. NELEM must be positive.");
    }
  else if (NERROR == 1107)
    {
      strcpy (REASON,
	      "Fatal Error in PEMATW. Abnormal termination of file ''pdcompos.pen''.");
    }
  else if (NERROR == 1108)
    {
      strcpy (REASON, "Fatal Error in PEMATW. Wrong atomic number.");
    }
  else if (NERROR == 1109)
    {
      strcpy (REASON,
	      "Fatal Error in PEMATW. STF (atoms/molecule) must be positive.");
    }
  else if (NERROR == 1110)
    {
      strcpy (REASON,
	      "Fatal Error in PEMATW. Element has been declared twice.");
    }
  else if (NERROR == 1111)
    {
      strcpy (REASON, "Fatal Error in PEMATW. The density must be positive.");
    }
  else if (NERROR == 1112)
    {
      strcpy (REASON, "Fatal Error in PEMATW. Too many shells.");
    }
  else if (NERROR == 1113)
    {
      strcpy (REASON, "Fatal Error in PEMATW. Wrong shell number.");
    }
  else if (NERROR == 1114)
    {
      strcpy (REASON, "Fatal Error in PEMATW. Unbalanced charges (element).");
    }
  else if (NERROR == 1115)
    {
      strcpy (REASON,
	      "Fatal Error in PEMATW. Too many oscillators. The parameter NOM should be increased.");
    }
  else if (NERROR == 1116)
    {
      strcpy (REASON,
	      "Fatal Error in PEMATW. Unbalanced charges (compound).");
    }
  else if (NERROR == 1117)
    {
      strcpy (REASON, "Fatal Error in PEMATW. FP is too large.");
    }
  else if (NERROR == 1118)
    {
      strcpy (REASON,
	      "Fatal Error in PEMATW. Inconsistent oscillator strength data.");
    }
  else if (NERROR == 1119)
    {
      strcpy (REASON,
	      "Fatal Error in PEMATW. Inconsistent oscillator-strength data.");
    }
  else if (NERROR == 1120)
    {
      strcpy (REASON,
	      "Fatal Error in PEMATW. Inconsistent oscillator-strength data (2).");
    }
  else if (NERROR == 1121)
    {
      strcpy (REASON,
	      "Fatal Error in PEMATW. Error in grouping the Compton profiles.");
    }
  else if (NERROR == 1201)
    {
      strcpy (REASON, "Fatal Error in START. E out of range.");
    }
  else if (NERROR == 1202)
    {
      strcpy (REASON, "Fatal Error in KNOCK: Incorrect particle type.");
    }
  else if (NERROR == 1203)
    {
      strcpy (REASON,
	      "Fatal Error in STORES: Not enough storage for secondary positrons.");
    }
  else if (NERROR == 1204)
    {
      strcpy (REASON, "Fatal Error in KNOCKF: Incorrect particle type.");
    }
  else if (NERROR == 1205)
    {
      strcpy (REASON, "Fatal Error in JUMPW: The particle is not a photon.");
    }
  else if (NERROR == 1206)
    {
      strcpy (REASON, "Fatal Error in KNOCKW: The particle is not a photon.");
    }
  else if (NERROR == 1301)
    {
      strcpy (REASON, "Fatal Error in EELa0. Negative total cross section.");
    }
  else if (NERROR == 1302)
    {
      strcpy (REASON, "Fatal Error in EELa0. Inconsistent arguments.");
    }
  else if (NERROR == 1303)
    {
      strcpy (REASON, "Fatal Error in EELaW. Wrong file.");
    }
  else if (NERROR == 1304)
    {
      strcpy (REASON, "Fatal Error in ESIaR. Corrupt material data file.");
    }
  else if (NERROR == 1305)
    {
      strcpy (REASON, "Fatal Error in ESIaR. Too many data points.");
    }
  else if (NERROR == 1306)
    {
      strcpy (REASON, "Fatal Error in ESIaR. Too many shells.");
    }
  else if (NERROR == 1307)
    {
      strcpy (REASON, "Fatal Error in ESIaR. Insufficient memory storage.");
    }
  else if (NERROR == 1308)
    {
      strcpy (REASON, "Fatal Error in ESIaW. Corrupt data file.");
    }
  else if (NERROR == 1309)
    {
      strcpy (REASON, "Fatal Error in ESIaW. Too many shells.");
    }
  else if (NERROR == 1310)
    {
      strcpy (REASON, "Fatal Error in PSIaR. Corrupt material data file.");
    }
  else if (NERROR == 1311)
    {
      strcpy (REASON, "Fatal Error in PSIaR. Too many data points.");
    }
  else if (NERROR == 1312)
    {
      strcpy (REASON, "Fatal Error in PSIaR. Too many shells.");
    }
  else if (NERROR == 1313)
    {
      strcpy (REASON, "Fatal Error in PSIaR. Insufficient memory storage.");
    }
  else if (NERROR == 1314)
    {
      strcpy (REASON, "Fatal Error in PSIaW. Corrupt data file.");
    }
  else if (NERROR == 1315)
    {
      strcpy (REASON, "Fatal Error in PSIaW. Too many shells.");
    }
  else if (NERROR == 1316)
    {
      strcpy (REASON, "Fatal Error in EBRaR. EBRR. Inconsistent format.");
    }
  else if (NERROR == 1317)
    {
      strcpy (REASON, "Fatal Error in EBRW. Corrupt file.");
    }
  else if (NERROR == 1318)
    {
      strcpy (REASON, "Fatal Error in EBRW. Check the bremss database file.");
    }
  else if (NERROR == 1319)
    {
      strcpy (REASON, "Fatal Error in RLMOM. Error code 0.");
    }
  else if (NERROR == 1320)
    {
      strcpy (REASON, "Fatal Error in RLMOM. Error code 1.");
    }
  else if (NERROR == 1321)
    {
      strcpy (REASON, "Fatal Error in RLMOM. Error code 2.");
    }
  else if (NERROR == 1322)
    {
      strcpy (REASON, "Fatal Error in RLMOM. Error code 3.");
    }
  else if (NERROR == 1323)
    {
      strcpy (REASON, "Fatal Error in RLMOM. Error code 4.");
    }
  else if (NERROR == 1324)
    {
      strcpy (REASON, "Fatal Error in BRAR. Inconsistent data.");
    }
  else if (NERROR == 1325)
    {
      strcpy (REASON,
	      "Fatal Error in BRAR. Corrupt data file (pdbrang.p08).");
    }
  else if (NERROR == 1326)
    {
      strcpy (REASON,
	      "Fatal Error in BRAW. Corrupt data file (pdbrang.p08).");
    }
  else if (NERROR == 1327)
    {
      strcpy (REASON, "Fatal Error in GRAaR. RITA initialisation error.");
    }
  else if (NERROR == 1328)
    {
      strcpy (REASON, "Fatal Error in GRAaR. RITA initialisation error.");
    }
  else if (NERROR == 1329)
    {
      strcpy (REASON,
	      "Fatal Error in GRAaR. RITA interpolation error is too large.");
    }
  else if (NERROR == 1330)
    {
      strcpy (REASON, "Fatal Error in GRAaW. Corrupt file.");
    }
  else if (NERROR == 1331)
    {
      strcpy (REASON, "Fatal Error in GRAaW. Corrupt file.");
    }
  else if (NERROR == 1332)
    {
      strcpy (REASON, "Fatal Error in GRAaW. Corrupt file.");
    }
  else if (NERROR == 1333)
    {
      strcpy (REASON, "Fatal Error in GPHaR. Corrupt material data file.");
    }
  else if (NERROR == 1334)
    {
      strcpy (REASON, "Fatal Error in GPHaR. Too many data points.");
    }
  else if (NERROR == 1335)
    {
      strcpy (REASON, "Fatal Error in GPHaR. Too many shells.");
    }
  else if (NERROR == 1336)
    {
      strcpy (REASON, "Fatal Error in GPHaR. Insufficient memory storage.");
    }
  else if (NERROR == 1337)
    {
      strcpy (REASON, "Fatal Error in GPHaW. Too many energies.");
    }
  else if (NERROR == 1338)
    {
      strcpy (REASON, "Fatal Error in GPHaW. Corrupt file.");
    }
  else if (NERROR == 1339)
    {
      strcpy (REASON, "Fatal Error in GPHaW. Too many shells.");
    }
  else if (NERROR == 1340)
    {
      strcpy (REASON, "Fatal Error in GPPaW. Corrupt file.");
    }
  else if (NERROR == 1341)
    {
      strcpy (REASON, "Fatal Error in RELAXR. NTRAN needs to be increased.");
    }
  else if (NERROR == 1342)
    {
      strcpy (REASON, "Fatal Error in RELAXR. Insufficient memory storage.");
    }
  else if (NERROR == 1343)
    {
      strcpy (REASON,
	      "Fatal Error in RELAXR. Negative transition probability?");
    }
  else if (NERROR == 1344)
    {
      strcpy (REASON, "Fatal Error in RELAXR. Rounding error is too large.");
    }
  else if (NERROR == 1345)
    {
      strcpy (REASON, "Fatal Error in RELAXW. The element is not loaded.");
    }
  else if (NERROR == 1346)
    {
      strcpy (REASON, "Fatal Error in RELAXW. NM needs to be increased.");
    }
  else if (NERROR == 1347)
    {
      strcpy (REASON,
	      "Fatal Error in EELdW. Electron cross section data are corrupt.");
    }
  else if (NERROR == 1348)
    {
      strcpy (REASON,
	      "Fatal Error in EELdW. Positron cross section data are corrupt.");
    }
  else if (NERROR == 1349)
    {
      strcpy (REASON,
	      "Fatal Error in EELdR. Error reading electron elastic DCS data.");
    }
  else if (NERROR == 1350)
    {
      strcpy (REASON,
	      "Fatal Error in EELdR. Electron cross section data are corrupt.");
    }
  else if (NERROR == 1351)
    {
      strcpy (REASON,
	      "Fatal Error in EELdR. Error reading positron elastic DCS data.");
    }
  else if (NERROR == 1352)
    {
      strcpy (REASON,
	      "Fatal Error in EELdR. Positron cross section data are corrupt.");
    }
  else if (NERROR == 1353)
    {
      strcpy (REASON, "Fatal Error in EELdR. RITA initialisation error.");
    }
  else if (NERROR == 1354)
    {
      strcpy (REASON, "Fatal Error in EELdR. RITA initialisation error.");
    }
  else if (NERROR == 1355)
    {
      strcpy (REASON, "Fatal Error in ELINIT. Corrupt data file.");
    }
  else if (NERROR == 1356)
    {
      strcpy (REASON, "Fatal Error in ELINIT. Corrupt data file.");
    }
  else if (NERROR == 1357)
    {
      strcpy (REASON, "Fatal Error in DCSEL. Energy out of range.");
    }
  else if (NERROR == 1401)
    {
      strcpy (REASON, "Fatal Error in IRND0. Negative point probability.");
    }
  else if (NERROR == 1402)
    {
      strcpy (REASON, "Fatal Error in RITAI0: N must be larger than 8.");
    }
  else if (NERROR == 1403)
    {
      strcpy (REASON, "Fatal Error in RITAI0: N must be less than NM=512.");
    }
  else if (NERROR == 1404)
    {
      strcpy (REASON,
	      "Fatal Error in RITAI0: XLOW must be larger than XHIGH.");
    }
  else if (NERROR == 1405)
    {
      strcpy (REASON, "Fatal Error in RITAI0: XLOW and NU are negative.");
    }
  else if (NERROR == 1901)
    {
      strcpy (REASON,
	      "Fatal Error in MERGE2. Increase the value of the parameter NP.");
    }
  else if (NERROR == 1902)
    {
      strcpy (REASON,
	      "Fatal Error in MERGE2. Increase the value of the parameter NP.");
    }
  else if (NERROR == 1903)
    {
      strcpy (REASON,
	      "Fatal Error in SORT2. Increase the value of the parameter NP.");
    }
  else if (NERROR == 1904)
    {
      strcpy (REASON, "Fatal Error in SPLINE: N is less than 4.");
    }
  else if (NERROR == 1905)
    {
      strcpy (REASON,
	      "Fatal Error in SPLINE: X values not in increasing order.");
    }
  else if (NERROR == 1906)
    {
      strcpy (REASON, "Fatal Error in SINTEG. Integral limits out of range.");
    }
  else if (NERROR == 1907)
    {
      strcpy (REASON, "Fatal Error int SLAG6: too few data points.");
    }
  else if (NERROR == 1908)
    {
      strcpy (REASON, "Fatal Error in RMOMX. Error code 1.");
    }
  else if (NERROR == 1909)
    {
      strcpy (REASON, "Fatal Error in RMOMX. Error code 2.");
    }
  else if (NERROR == 1910)
    {
      strcpy (REASON, "Fatal Error in RMOMX. Error code 3.");
    }
  else if (NERROR == 1911)
    {
      strcpy (REASON, "Fatal Error in RMOMX. Error code 4.");
    }
  else if (NERROR == 1912)
    {
      strcpy (REASON, "Fatal Error in RNDG30: initialisation error (1).");
    }
  else if (NERROR == 1913)
    {
      strcpy (REASON, "Fatal Error in RNDG30: initialisation error (2).");
    }
  else if (NERROR == 2001)
    {
      strcpy (REASON,
	      "Fatal Error in GEOMIN: The input and output units must be different.");
    }
  else if (NERROR == 2002)
    {
      strcpy (REASON, "Fatal Error in GEOMIN: What do you mean?");
    }
  else if (NERROR == 2003)
    {
      strcpy (REASON, "Fatal Error in GEOMIN: Incorrect label format.");
    }
  else if (NERROR == 2004)
    {
      strcpy (REASON, "Fatal Error in GEOMIN: Same label for two surfaces.");
    }
  else if (NERROR == 2005)
    {
      strcpy (REASON,
	      "Fatal Error in GEOMIN: The parameter NS must be increased.");
    }
  else if (NERROR == 2006)
    {
      strcpy (REASON,
	      "Fatal Error in GEOMIN: Incorrect format of surface indices.");
    }
  else if (NERROR == 2007)
    {
      strcpy (REASON, "Fatal Error in GEOMIN: Incorrect surface indices.");
    }
  else if (NERROR == 2008)
    {
      strcpy (REASON,
	      "Fatal Error in GEOMIN: NPINP is too small (check PARINP).");
    }
  else if (NERROR == 2009)
    {
      strcpy (REASON,
	      "Fatal Error in GEOMIN: Scale factor less than 1.0E-15.)");
    }
  else if (NERROR == 2010)
    {
      strcpy (REASON,
	      "Fatal Error in GEOMIN: Scale factor less than 1.0E-15.)");
    }
  else if (NERROR == 2011)
    {
      strcpy (REASON,
	      "Fatal Error in GEOMIN: Scale factor less than 1.0E-15.)");
    }
  else if (NERROR == 2012)
    {
      strcpy (REASON, "Fatal Error in GEOMIN: What do you mean?");
    }
  else if (NERROR == 2013)
    {
      strcpy (REASON,
	      "Fatal Error in GEOMIN: NPINP is too small (check PARINP)");
    }
  else if (NERROR == 2014)
    {
      strcpy (REASON, "Fatal Error in GEOMIN: What do you mean?");
    }
  else if (NERROR == 2015)
    {
      strcpy (REASON,
	      "Fatal Error in GEOMIN: NPINP is too small (check PARINP)");
    }
  else if (NERROR == 2016)
    {
      strcpy (REASON, "Fatal Error in GEOMIN: What do you mean?");
    }
  else if (NERROR == 2017)
    {
      strcpy (REASON, "Fatal Error in GEOMIN: Incorrect label format.");
    }
  else if (NERROR == 2018)
    {
      strcpy (REASON,
	      "Fatal Error in GEOMIN: Same label for two bodies (or modules).");
    }
  else if (NERROR == 2019)
    {
      strcpy (REASON,
	      "Fatal Error in GEOMIN: The parameter NB must be increased.");
    }
  else if (NERROR == 2020)
    {
      strcpy (REASON,
	      "Fatal Error in GEOMIN: Incorrect material definition line.");
    }
  else if (NERROR == 2021)
    {
      strcpy (REASON, "Fatal Error in GEOMIN: Undefined surface label.");
    }
  else if (NERROR == 2022)
    {
      strcpy (REASON,
	      "Fatal Error in GEOMIN: The last limiting surface has been defined twice.");
    }
  else if (NERROR == 2023)
    {
      strcpy (REASON,
	      "Fatal Error in GEOMIN: The number of limiting surfaces is too large.");
    }
  else if (NERROR == 2024)
    {
      strcpy (REASON, "Fatal Error in GEOMIN: Check side pointer value.");
    }
  else if (NERROR == 2025)
    {
      strcpy (REASON, "Fatal Error in GEOMIN: Undefined body label.");
    }
  else if (NERROR == 2026)
    {
      strcpy (REASON, "Fatal Error in GEOMIN: This body is a module.");
    }
  else if (NERROR == 2027)
    {
      strcpy (REASON,
	      "Fatal Error in GEOMIN: The number of limiting surfaces is too large.");
    }
  else if (NERROR == 2028)
    {
      strcpy (REASON, "Fatal Error in GEOMIN: Undefined body label.");
    }
  else if (NERROR == 2029)
    {
      strcpy (REASON, "Fatal Error in GEOMIN: This module is a body.");
    }
  else if (NERROR == 2030)
    {
      strcpy (REASON,
	      "Fatal Error in GEOMIN: The number of limiting surfaces is too large.");
    }
  else if (NERROR == 2031)
    {
      strcpy (REASON, "Fatal Error in GEOMIN: What do you mean?");
    }
  else if (NERROR == 2032)
    {
      strcpy (REASON, "Fatal Error in GEOMIN: Incorrect label format.");
    }
  else if (NERROR == 2033)
    {
      strcpy (REASON,
	      "Fatal Error in GEOMIN: Same label for two bodies (or modules).");
    }
  else if (NERROR == 2034)
    {
      strcpy (REASON,
	      "Fatal Error in GEOMIN: The parameter NB must be increased.");
    }
  else if (NERROR == 2035)
    {
      strcpy (REASON,
	      "Fatal Error in GEOMIN: Incorrect material definition line.");
    }
  else if (NERROR == 2036)
    {
      strcpy (REASON, "Fatal Error in GEOMIN: Undefined surface label.");
    }
  else if (NERROR == 2037)
    {
      strcpy (REASON,
	      "Fatal Error in GEOMIN: The last limiting surface has been defined twice.");
    }
  else if (NERROR == 2038)
    {
      strcpy (REASON,
	      "Fatal Error in GEOMIN: The number of limiting surfaces is too large.");
    }
  else if (NERROR == 2039)
    {
      strcpy (REASON, "Fatal Error in GEOMIN: Check side pointer value.");
    }
  else if (NERROR == 2040)
    {
      strcpy (REASON, "Fatal Error in GEOMIN: Undefined body label.");
    }
  else if (NERROR == 2041)
    {
      strcpy (REASON, "Fatal Error in GEOMIN: This body is a module.");
    }
  else if (NERROR == 2042)
    {
      strcpy (REASON,
	      "Fatal Error in GEOMIN: You are trying to assign two mothers to the last body.");
    }
  else if (NERROR == 2043)
    {
      strcpy (REASON,
	      "Fatal Error in GEOMIN: The number of limiting surfaces is too large.");
    }
  else if (NERROR == 2044)
    {
      strcpy (REASON, "Fatal Error in GEOMIN: Undefined body label.");
    }
  else if (NERROR == 2045)
    {
      strcpy (REASON, "Fatal Error in GEOMIN: This module is a body.");
    }
  else if (NERROR == 2046)
    {
      strcpy (REASON,
	      "Fatal Error in GEOMIN: You are trying to assign two mothers to the last module.");
    }
  else if (NERROR == 2047)
    {
      strcpy (REASON,
	      "Fatal Error in GEOMIN: The number of limiting surfaces is too large.");
    }
  else if (NERROR == 2048)
    {
      strcpy (REASON, "Fatal Error in GEOMIN: What do you mean?");
    }
  else if (NERROR == 2049)
    {
      strcpy (REASON,
	      "Fatal Error in GEOMIN: NPINP is too small (check PARINP)");
    }
  else if (NERROR == 2050)
    {
      strcpy (REASON, "Fatal Error in GEOMIN: What do you mean?");
    }
  else if (NERROR == 2051)
    {
      strcpy (REASON, "Fatal Error in GEOMIN: Incorrect label format.");
    }
  else if (NERROR == 2052)
    {
      strcpy (REASON,
	      "Fatal Error in GEOMIN: Same label for two bodies or modules.");
    }
  else if (NERROR == 2053)
    {
      strcpy (REASON,
	      "Fatal Error in GEOMIN: The cloned object must be a module.");
    }
  else if (NERROR == 2054)
    {
      strcpy (REASON, "Fatal Error in GEOMIN: This module is not defined.");
    }
  else if (NERROR == 2055)
    {
      strcpy (REASON,
	      "Fatal Error in GEOMIN: The selected object is a body.");
    }
  else if (NERROR == 2056)
    {
      strcpy (REASON,
	      "Fatal Error in GEOMIN: The label does not correspond to a module.");
    }
  else if (NERROR == 2057)
    {
      strcpy (REASON, "Fatal Error in GEOMIN: What do you mean?");
    }
  else if (NERROR == 2058)
    {
      strcpy (REASON,
	      "Fatal Error in GEOMIN: NPINP is too small (check PARINP)");
    }
  else if (NERROR == 2059)
    {
      strcpy (REASON, "Fatal Error in GEOMIN: What do you mean?");
    }
  else if (NERROR == 2060)
    {
      strcpy (REASON,
	      "Fatal Error in GEOMIN: The parameter NS must be increased.");
    }
  else if (NERROR == 2061)
    {
      strcpy (REASON,
	      "Fatal Error in GEOMIN: The parameter NB must be increased.");
    }
  else if (NERROR == 2062)
    {
      strcpy (REASON, "Fatal Error in GEOMIN: Something wrong...");
    }
  else if (NERROR == 2063)
    {
      strcpy (REASON,
	      "Fatal Error in GEOMIN: The limiting body or module is not yet defined");
    }
  else if (NERROR == 2064)
    {
      strcpy (REASON,
	      "Fatal Error in GEOMIN: The limiting body or module is not yet defined");
    }
  else if (NERROR == 2065)
    {
      strcpy (REASON, "Fatal Error in GEOMIN: Incorrect label format.");
    }
  else if (NERROR == 2066)
    {
      strcpy (REASON, "Fatal Error in GEOMIN: What do you mean?");
    }
  else if (NERROR == 2067)
    {
      strcpy (REASON, "Fatal Error in GEOMIN: Too many include levels.");
    }
  else if (NERROR == 2068)
    {
      strcpy (REASON, "Fatal Error in GEOMIN: Too many include levels.");
    }
  else if (NERROR == 2069)
    {
      strcpy (REASON,
	      "Fatal Error in GEOMIN: The parameter NB must be increased.");
    }
  else if (NERROR == 2070)
    {
      strcpy (REASON,
	      "Fatal Error in GEOMIN: The parameter NS must be increased.");
    }
  else if (NERROR == 2071)
    {
      strcpy (REASON,
	      "Fatal Error in GEOMIN: The parameter NXG is too small.");
    }
  else if (NERROR == 2072)
    {
      strcpy (REASON,
	      "Fatal Error in GEOMIN: Possibly unresolved body or module.");
    }
  else if (NERROR == 2073)
    {
      strcpy (REASON, "Fatal Error in GEOMIN: Inconsistent body label.");
    }
  else if (NERROR == 2074)
    {
      strcpy (REASON, "Fatal Error in GEOMIN: Inconsistent side pointers.");
    }
  else if (NERROR == 2075)
    {
      strcpy (REASON, "Fatal Error in GEOMIN: Wrong input format.");
    }
  //ERRORS PENEASY
  else if (NERROR == 146)
    {
      strcpy (REASON,
	      "samplePosition:ERROR: sampling source efficiency is lower than 0.1%%.");
    }
  else if (NERROR == 147)
    {
      strcpy (REASON, "iniconfig:ERROR: incorrect section header;");
    }
  else if (NERROR == 148)
    {
      strcpy (REASON, "iniconfig:ERROR: too many requested histories");
    }
  else if (NERROR == 149)
    {
      strcpy (REASON, "iniconfig:ERROR: refresh interval must be positive.");
    }
  else if (NERROR == 150)
    {
      strcpy (REASON, "iniconfig:ERROR: Invalid RNG seeds.");
    }
  else if (NERROR == 151)
    {
      strcpy (REASON, "iniconfig:ERROR: unable to open seeds file.");
    }
  else if (NERROR == 152)
    {
      strcpy (REASON, "iniconfig:ERROR: unable to open restart file.");
    }
  else if (NERROR == 153)
    {
      strcpy (REASON,
	      "iniconfig:ERROR: interval between dumps must be positive.");
    }
  else if (NERROR == 154)
    {
      strcpy (REASON, "iniconfig:ERROR: End-Of-Section mark not found");
    }
  else if (NERROR == 155)
    {
      strcpy (REASON, "inigeo:ERROR: incorrect section header");
    }
  else if (NERROR == 156)
    {
      strcpy (REASON, "inigeo:ERROR: unable to open quadrics file");
    }
  else if (NERROR == 157)
    {
      strcpy (REASON,
	      "inigeo:ERROR: too many materials; enlarge MAXMAT parameter and recompile.");
    }
  else if (NERROR == 158)
    {
      strcpy (REASON, "inigeo:ERROR: no geometry defined.");
    }
  else if (NERROR == 159)
    {
      strcpy (REASON, "inigeo:ERROR: End-Of-Section mark not found");
    }
  else if (NERROR == 160)
    {
      strcpy (REASON, "inipen:ERROR: incorrect section header;");
    }
  else if (NERROR == 161)
    {
      strcpy (REASON, "inipen:ERROR: Max number of materials exceeded");
    }
  else if (NERROR == 162)
    {
      strcpy (REASON, "inipen:ERROR: Invalid MAT index");
    }
  else if (NERROR == 163)
    {
      strcpy (REASON, "inipen:ERROR: incomplete list of parameters for MAT");
    }
  else if (NERROR == 164)
    {
      strcpy (REASON,
	      "inipen:ERROR: DSMAX must be larger than zero even if electrons are not transported.");
    }
  else if (NERROR == 165)
    {
      strcpy (REASON,
	      "inipen:ERROR: There are more materials declared in the geometry file than defined in the config file.");
    }
  else if (NERROR == 166)
    {
      strcpy (REASON, "inipen:ERROR: End-Of-Section mark not found");
    }
  else if (NERROR == 167)
    {
      strcpy (REASON, "getline:ERROR: unable to read line.");
    }
  else if (NERROR == 168)
    {
      strcpy (REASON, "seeki:ERROR: value outside range, xc>x(n):");
    }
  else if (NERROR == 169)
    {
      strcpy (REASON, "seeki:ERROR: value outside range, xc<x(1):");
    }
  else if (NERROR == 170)
    {
      strcpy (REASON, "BIGSinisrc:ERROR: incorrect section header;");
    }
  else if (NERROR == 171)
    {
      strcpy (REASON, "BIGSinisrc:ERROR: Unable to find End-Of-Section");
    }
  else if (NERROR == 172)
    {
      strcpy (REASON, "BIGSinisrc:ERROR: expecting to find ON or OFF");
    }
  else if (NERROR == 173)
    {
      strcpy (REASON, "BIGSinisrc:ERROR: invalid particle type");
    }
  else if (NERROR == 174)
    {
      strcpy (REASON,
	      "BIGSinisrc:ERROR: invalid polarization switch, should be 0 or 1");
    }
  else if (NERROR == 175)
    {
      strcpy (REASON,
	      "BIGSinisrc:ERROR: vector P={P1,P2,P3} must be P^2 <= 1");
    }
  else if (NERROR == 176)
    {
      strcpy (REASON, "BIGSinisrc:ERROR: null direction.");
    }
  else if (NERROR == 177)
    {
      strcpy (REASON, "BIGSinisrc:ERROR: theta1 is less than theta0.");
    }
  else if (NERROR == 178)
    {
      strcpy (REASON, "BIGSinisrc:ERROR: invalid interval.");
    }
  else if (NERROR == 179)
    {
      strcpy (REASON, "BIGSinisrc:ERROR: Invalid entry. Must be 0 or 1.");
    }
  else if (NERROR == 180)
    {
      strcpy (REASON, "BIGSinisrc:ERROR: unable to open spectrum file.");
    }
  else if (NERROR == 181)
    {
      strcpy (REASON, "BIGSinisrc:ERROR: invalid entry");
    }
  else if (NERROR == 182)
    {
      strcpy (REASON, "BIGSinisrc:ERROR: negative energy");
    }
  else if (NERROR == 183)
    {
      strcpy (REASON, "BIGSinisrc:ERROR: decreasing energy");
    }
  else if (NERROR == 184)
    {
      strcpy (REASON, "BIGSinisrc:ERROR: too many bins in spectrum;");
      strcpy (REASON, "              enlarge NEMAX");
    }
  else if (NERROR == 185)
    {
      strcpy (REASON, "BIGSinisrc:ERROR: at least 1 bin must be defined");
    }
  else if (NERROR == 186)
    {
      strcpy (REASON, "BIGSinisrc:ERROR: all probabilities are zero.");
    }
  else if (NERROR == 187)
    {
      strcpy (REASON, "BIGSinisrc:ERROR: End-Of-Section mark not found");
    }
  else if (NERROR == 188)
    {
      strcpy (REASON,
	      "inisource:ERROR: PSF source ON is incompatible with other sources ON.");
    }
  else if (NERROR == 189)
    {
      strcpy (REASON, "EDPinitally:ERROR: incorrect section header;");
    }
  else if (NERROR == 190)
    {
      strcpy (REASON,
	      "EDPinitally:ERROR: Unable to find End-Of-Section mark");
    }
  else if (NERROR == 191)
    {
      strcpy (REASON, "EDPinitally:ERROR: expecting to find ON or OFF");
    }
  else if (NERROR == 192)
    {
      strcpy (REASON, "iniforce:ERROR: incorrect section header");
    }
  else if (NERROR == 193)
    {
      strcpy (REASON, "iniforce:ERROR: Unable to find End-Of-Section mark");
    }
  else if (NERROR == 194)
    {
      strcpy (REASON, "iniforce:ERROR: expecting to find ON or OFF");
    }
  else if (NERROR == 195)
    {
      strcpy (REASON,
	      "iniforce:ERROR: unable to read line containing: MAT,KPAR,ICOL,forcing");
    }
  else if (NERROR == 196)
    {
      strcpy (REASON, "iniforce:ERROR: invalid MAT");
    }
  else if (NERROR == 197)
    {
      strcpy (REASON, "iniforce:ERROR: proton forcing not implemented yet");
    }
  else if (NERROR == 198)
    {
      strcpy (REASON, "iniforce:ERROR: KPAR must be in [1,3]");
    }
  else if (NERROR == 199)
    {
      strcpy (REASON, "iniforce:ERROR: ICOL must be in [0,8]");
    }
  else if (NERROR == 200)
    {
      strcpy (REASON, "iniforce:ERROR: FORCING must not be < 1");
    }
  else if (NERROR == 201)
    {
      strcpy (REASON, "iniforce:ERROR: End-Of-Section mark not found");
    }
  else if (NERROR == 202)
    {
      strcpy (REASON, "inisplit:ERROR: incorrect section header");
    }
  else if (NERROR == 203)
    {
      strcpy (REASON, "inisplit:ERROR: Unable to find End-Of-Section mark");
    }
  else if (NERROR == 204)
    {
      strcpy (REASON, "inisplit:ERROR: expecting to find ON or OFF");
    }
  else if (NERROR == 205)
    {
      strcpy (REASON, "inisplit:ERROR: Min weight must be >0.");
    }
  else if (NERROR == 206)
    {
      strcpy (REASON,
	      "inisplit:ERROR: Photon polarization is active. Only simple splitting can be used");
    }
  else if (NERROR == 207)
    {
      strcpy (REASON, "inisplit:ERROR: Invalid factor");
    }
  else if (NERROR == 208)
    {
      strcpy (REASON, "inisplit:ERROR: Invalid sign");
    }
  else if (NERROR == 209)
    {
      strcpy (REASON, "inisplit:ERROR: Invalid sector");
    }
  else if (NERROR == 210)
    {
      strcpy (REASON, "inisplit:ERROR: End-Of-Section mark not found");
    }
  else if (NERROR == 211)
    {
      strcpy (REASON, "splitting:ERROR: internal error");
    }
  else if (NERROR == 212)
    {
      strcpy (REASON, "inirussia:ERROR: incorrect section header");
    }
  else if (NERROR == 213)
    {
      strcpy (REASON, "inirussia:ERROR: Unable to find End-Of-Section mark");
    }
  else if (NERROR == 214)
    {
      strcpy (REASON, "inirussia:ERROR: expecting to find ON or OFF");
    }
  else if (NERROR == 215)
    {
      strcpy (REASON, "inirussia:ERROR: Invalid survival probability");
    }
  else if (NERROR == 216)
    {
      strcpy (REASON, "inirussia:ERROR: End-Of-Section mark not found");
    }
  else if (NERROR == 217)
    {
      strcpy (REASON, "SDDinitially:ERROR: incorrect section header");
    }
  else if (NERROR == 218)
    {
      strcpy (REASON, "inivox:ERROR: invalid transparent mat.");
    }
  else if (NERROR == 219)
    {
      strcpy (REASON,
	      "inivox:ERROR: there must be ONE body with transparent material.");
    }
  else if (NERROR == 220)
    {
      strcpy (REASON,
	      "inivox:ERROR: granularity must be between 2 and maxGranul");
    }
  else if (NERROR == 221)
    {
      strcpy (REASON,
	      "inivox:ERROR: internal inconsistency; expecting vacuum");
    }
  else if (NERROR == 222)
    {
      strcpy (REASON, "readvox:ERROR: unable to open voxels file");
    }
  else if (NERROR == 223)
    {
      strcpy (REASON, "getvoxline:ERROR: unable to read vox file");
    }
  else if (NERROR == 224)
    {
      strcpy (REASON, "readvox:ERROR: incorrect section header");
    }
  else if (NERROR == 225)
    {
      strcpy (REASON, "readvox:ERROR: invalid no. voxels.");
    }
  else if (NERROR == 226)
    {
      strcpy (REASON, "readvox:ERROR: No. voxels exceeds nvoxmax.");
    }
  else if (NERROR == 227)
    {
      strcpy (REASON, "readvox:ERROR: not enough memory.");
    }
  else if (NERROR == 228)
    {
      strcpy (REASON, "inimassvox:ERROR: not enough memory for temp arrays.");
    }
  else if (NERROR == 229)
    {
      strcpy (REASON, "inimassvox:ERROR: not enough memory.");
    }
  else if (NERROR == 230)
    {
      strcpy (REASON, "readvox:ERROR: voxel side too small.");
    }
  else if (NERROR == 231)
    {
      strcpy (REASON, "readvox:ERROR: VBB too large.");
    }
  else if (NERROR == 232)
    {
      strcpy (REASON,
	      "readvox:ERROR: column numbers must be between 1 and maxCol.");
    }
  else if (NERROR == 233)
    {
      strcpy (REASON, "readvox:ERROR: End-Of-Section mark not found.");
    }
  else if (NERROR == 234)
    {
      strcpy (REASON, "readvox:ERROR: invalid entry at line line.");
    }
  else if (NERROR == 235)
    {
      strcpy (REASON, "readvox:ERROR: Line should be blank, line line.");
    }
  else if (NERROR == 236)
    {
      strcpy (REASON, "writeMassvox:ERROR: unable to open file to write.");
    }
  else if (NERROR == 237)
    {
      strcpy (REASON,
	      "SDDinitally:ERROR: Unable to find End-Of-Section mark.");
    }
  else if (NERROR == 238)
    {
      strcpy (REASON, "SDDinitally:ERROR: expecting to find ON or OFF.");
    }
  else if (NERROR == 239)
    {
      strcpy (REASON, "SDDinitally:ERROR: invalid entry.");
    }
  else if (NERROR == 240)
    {
      strcpy (REASON, "SDDinitally:ERROR: Zero bins defined.");
    }
  else if (NERROR == 241)
    {
      strcpy (REASON,
	      "SDDinitally:ERROR: nbinmax max no. of megabins exceeded.");
    }
  else if (NERROR == 242)
    {
      strcpy (REASON, "SDDinitally:ERROR: expecting 1,0 or -1.");
    }
  else if (NERROR == 243)
    {
      strcpy (REASON, "SDDinitally:ERROR: not enough memory.");
    }
  else if (NERROR == 244)
    {
      strcpy (REASON, "SDDinitally:ERROR: End-Of-Section mark not found.");
    }
  else if (NERROR == 245)
    {
      strcpy (REASON, "PSFinitally:ERROR: incorrect section header.");
    }
  else if (NERROR == 246)
    {
      strcpy (REASON,
	      "PSFinitally:ERROR: Unable to find End-Of-Section mark.");
    }
  else if (NERROR == 247)
    {
      strcpy (REASON, "PSFinitally:ERROR: expecting to find ON or OFF.");
    }
  else if (NERROR == 248)
    {
      strcpy (REASON,
	      "PSFinitally:ERROR: IAEA PSF format requested but not available.");
    }
  else if (NERROR == 249)
    {
      strcpy (REASON, "PSFinitally:ERROR: PSF format must be 0 or 1.");
    }
  else if (NERROR == 250)
    {
      strcpy (REASON,
	      "PSFinitally:ERROR: detection material out of range: 1,MAXMAT.");
    }
  else if (NERROR == 251)
    {
      strcpy (REASON,
	      "PSFinitally:ERROR: PSF detection material must be a perfect absorbent.");
    }
  else if (NERROR == 252)
    {
      strcpy (REASON, "PSFinitally:ERROR: Could not append PSF.");
    }
  else if (NERROR == 253)
    {
      strcpy (REASON, "PSFinitally:ERROR: End-Of-Section mark not found.");
    }
  else if (NERROR == 254)
    {
      strcpy (REASON, "IAEAiniwrite:ERROR: Unable to open PSF.");
    }
  else if (NERROR == 255)
    {
      strcpy (REASON,
	      "IAEAiniwrite:ERROR: Could not open PSF to append; make sure the file exists and it is accessible.");
    }
  else if (NERROR == 256)
    {
      strcpy (REASON, "IAEAiniwrite:ERROR: Unable to set an extra variable.");
    }
  else if (NERROR == 257)
    {
      strcpy (REASON,
	      "IAEAiniwrite:internalERROR: Extra variable index is out of range.");
    }
  else if (NERROR == 258)
    {
      strcpy (REASON,
	      "IAEAiniwrite:internalERROR: Extra variable type is out of range.");
    }
  else if (NERROR == 259)
    {
      strcpy (REASON,
	      "IAEAiniwrite:internalERROR: Undefined error while setting an extra variable.");
    }
  else if (NERROR == 260)
    {
      strcpy (REASON, "IAEAwrite:ERROR: protons not implemented.");
    }
  else if (NERROR == 261)
    {
      strcpy (REASON, "IAEAwrite:ERROR: Invalid KPAR.");
    }
  else if (NERROR == 262)
    {
      strcpy (REASON, "IAEAwrite:ERROR: Unable to write particle.");
    }
  else if (NERROR == 263)
    {
      strcpy (REASON, "VDDinitially:ERROR: incorrect section header");
    }
  else if (NERROR == 264)
    {
      strcpy (REASON,
	      "VDDinitially:ERROR: voxel dose tally is ON but no voxelized geometry has been defined.");
    }
  else if (NERROR == 265)
    {
      strcpy (REASON,
	      "VDDinitally:ERROR: Unable to find End-Of-Section mark.");
    }
  else if (NERROR == 266)
    {
      strcpy (REASON, "VDDinitally:ERROR: expecting to find ON or OFF.");
    }
  else if (NERROR == 267)
    {
      strcpy (REASON, "VDDinitally:ERROR: invalid ROI.");
    }
  else if (NERROR == 268)
    {
      strcpy (REASON, "VDDinitally:ERROR: expecting 1,0 or -1.");
    }
  else if (NERROR == 269)
    {
      strcpy (REASON, "VDDinitally:ERROR: End-Of-Section mark not found.");
    }
  else if (NERROR == 270)
    {
      strcpy (REASON, "VDDinitally:ERROR: not enough memory.");
    }
  else if (NERROR == 271)
    {
      strcpy (REASON, "CDDinitally:ERROR: incorrect section header.");
    }
  else if (NERROR == 272)
    {
      strcpy (REASON,
	      "CDDinitally:ERROR: Unable to find End-Of-Section mark.");
    }
  else if (NERROR == 273)
    {
      strcpy (REASON, "CDDinitally:ERROR: expecting to find ON or OFF.");
    }
  else if (NERROR == 274)
    {
      strcpy (REASON, "CDDinitally:ERROR: invalid entry.");
    }
  else if (NERROR == 275)
    {
      strcpy (REASON, "CDDinitally:ERROR: Too many bins.");
    }
  else if (NERROR == 276)
    {
      strcpy (REASON, "CDDinitally:ERROR: End-Of-Section mark not found.");
    }
  else if (NERROR == 277)
    {
      strcpy (REASON, "PSFinisrc:ERROR: incorrect section header.");
    }
  else if (NERROR == 278)
    {
      strcpy (REASON, "PSFinisrc:ERROR: Unable to find End-Of-Section mark.");
    }
  else if (NERROR == 279)
    {
      strcpy (REASON, "PSFinisrc:ERROR: expecting to find ON or OFF.");
    }
  else if (NERROR == 280)
    {
      strcpy (REASON,
	      "PSFinisrc:ERROR: IAEA PSF format requested but not available.");
    }
  else if (NERROR == 281)
    {
      strcpy (REASON, "PSFinisrc:ERROR: PSF format must be 0 or 1.");
    }
  else if (NERROR == 282)
    {
      strcpy (REASON, "PSFinisrc:ERROR: split < 1.");
    }
  else if (NERROR == 283)
    {
      strcpy (REASON,
	      "PSFinisrc:ERROR: invalid entry. VALIDATE field must be 0 or 1.");
    }
  else if (NERROR == 284)
    {
      strcpy (REASON, "PSFinisrc:ERROR: cannot open the PSF.");
    }
  else if (NERROR == 285)
    {
      strcpy (REASON,
	      "PSFinisrc:ERROR: No. of histories in PSF exceeds nmax.");
    }
  else if (NERROR == 286)
    {
      strcpy (REASON, "PSFinisrc:ERROR: invalid KPAR: at line:.");
    }
  else if (NERROR == 287)
    {
      strcpy (REASON, "PSFinisrc:ERROR: invalid energy(eV): at line:.");
    }
  else if (NERROR == 288)
    {
      strcpy (REASON,
	      "PSFinisrc:ERROR: null vector direction found at line:.");
    }
  else if (NERROR == 289)
    {
      strcpy (REASON,
	      "PSFinisrc:ERROR: PSF end-of-file reached, not enough particles to initialize.");
    }
  else if (NERROR == 290)
    {
      strcpy (REASON,
	      "PSFinisrc:ERROR: Inconsistency found in PSF incremental history no.");
    }
  else if (NERROR == 291)
    {
      strcpy (REASON, "PSFinisrc:ERROR: End-Of-Section mark not found.");
    }
  else if (NERROR == 292)
    {
      strcpy (REASON, "getpar:ERROR: unable to read PSF line no.:.");
    }
  else if (NERROR == 293)
    {
      strcpy (REASON, "getpar:ERROR: invalid or missing datum in PSF line:.");
    }
  else if (NERROR == 294)
    {
      strcpy (REASON, "checkFormat:ERROR: cannot open the PSF.\n");
    }
  else if (NERROR == 295)
    {
      strcpy (REASON, "ckeckFormat:ERROR: unable to read first PSF line.\n");
    }
  else if (NERROR == 296)
    {
      strcpy (REASON, "checkFormat:ERROR: unable to identify PSF format.\n");
    }
  else if (NERROR == 297)
    {
      strcpy (REASON, "getparIAEA:ERROR: unable to read particle no.:.\n");
    }
  else if (NERROR == 298)
    {
      strcpy (REASON,
	      "getparIAEA:ERROR: EOF reached when attempting to read particle no.:.\n");
    }
  else if (NERROR == 299)
    {
      strcpy (REASON,
	      "getparIAEA:ERROR: undefined error while attempting to read particle no.:.\n");
    }
  else if (NERROR == 300)
    {
      strcpy (REASON,
	      "getparIAEA:ERROR: Invalid KPAR: found at particle no.:\n");
    }
  else if (NERROR == 301)
    {
      strcpy (REASON, "IAEAiniread:ERROR: unable to open PSF.\n");
    }
  else if (NERROR == 302)
    {
      strcpy (REASON,
	      "IAEAiniread:ERROR: Unable to get the total number of particles.\n");
    }
  else if (NERROR == 303)
    {
      strcpy (REASON,
	      "IAEAiniread:ERROR: Unable to get the total number of histories.\n");
    }
  else if (NERROR == 304)
    {
      strcpy (REASON,
	      "IAEAiniread:ERROR: No. of histories in PSF exceeds nmax.\n");
    }
  else if (NERROR == 305)
    {
      strcpy (REASON,
	      "IAEAiniread:ERROR: Unable to get constant variables.\n");
    }
  else if (NERROR == 306)
    {
      strcpy (REASON,
	      "IAEAiniread:ERROR: Variable constant index out of range.\n");
    }
  else if (NERROR == 307)
    {
      strcpy (REASON,
	      "IAEAiniread:ERROR: More than xtraMax extra int or real variables in PSF.\n");
    }
  else if (NERROR == 308)
    {
      strcpy (REASON,
	      "IAEAiniread:ERROR: Unable to check file size and byte order of the PSF.\n");
    }
  else if (NERROR == 309)
    {
      strcpy (REASON, "IAEAiniread:ERROR: The function fseek fails.\n");
    }
  else if (NERROR == 310)
    {
      strcpy (REASON,
	      "IAEAiniread:ERROR: File size inconsistent with header checksum.\n");
    }
  else if (NERROR == 311)
    {
      strcpy (REASON, "IAEAiniread:ERROR: There is a byte order mismatch.\n");
    }
  else if (NERROR == 312)
    {
      strcpy (REASON,
	      "IAEAiniread:ERROR: There is a file size and byte order mismatch.\n");
    }
  else if (NERROR == 313)
    {
      strcpy (REASON, "IAEAiniread:ERROR: Unidentified error code.\n");
    }
  else if (NERROR == 314)
    {
      strcpy (REASON, "IAEAiniread:internalERROR: invalid KPAR.\n");
    }
  else if (NERROR == 315)
    {
      strcpy (REASON,
	      "IAEAiniread:ERROR: energy out of range found in particle no.:.\n");
    }
  else if (NERROR == 316)
    {
      strcpy (REASON,
	      "IAEAiniread:ERROR: null direction vector found in particle no.:.\n");
    }
  else if (NERROR == 317)
    {
      strcpy (REASON,
	      "IAEAiniread:ERROR: No. of particles inferred from PSF differs from that stated in header file.\n");
    }
  else if (NERROR == 318)
    {
      strcpy (REASON,
	      "IAEAiniread:ERROR: No. of histories inferred from PSF is greater than that stated in header file.\n");
    }
  else if (NERROR == 319)
    {
      strcpy (REASON, "IAEAiniread:ERROR: Unable to get maximum energy.\n");
    }
  else if (NERROR == 320)
    {
      strcpy (REASON, "IAEAiniread:ERROR: unable to close PSF.\n");
    }
  else if (NERROR == 321)
    {
      strcpy (REASON, "IAEAiniread:ERROR: unable to open PSF.\n");
    }
  else if (NERROR == 322)
    {
      strcpy (REASON,
	      "IAEAiniread:ERROR: PSF end-of-file reached, not enough particles to initialize.\n");
    }
  else if (NERROR == 323)
    {
      strcpy (REASON,
	      "IAEAiniread:ERROR: Inconsistency found in PSF incremental history.\n");
    }
  else if (NERROR == 324)
    {
      strcpy (REASON, "PSFinitally:ERROR: incorrect section header.\n");
    }
  else if (NERROR == 325)
    {
      strcpy (REASON,
	      "PSFinitally:ERROR: Unable to find End-Of-Section mark.\n");
    }
  else if (NERROR == 326)
    {
      strcpy (REASON, "PSFinitally:ERROR: expecting to find ON or OFF.\n");
    }
  else if (NERROR == 327)
    {
      strcpy (REASON,
	      "PSFinitally:ERROR: IAEA PSF format requested but not available.\n");
    }
  else if (NERROR == 328)
    {
      strcpy (REASON, "PSFinitally:ERROR: PSF format must be 0 or 1.\n");
    }
  else if (NERROR == 329)
    {
      strcpy (REASON,
	      "PSFinitally:ERROR: detection material out of range.\n");
    }
  else if (NERROR == 330)
    {
      strcpy (REASON,
	      "PSFinitally:ERROR: PSF detection material must be a perfect absorbent.\n");
    }
  else if (NERROR == 331)
    {
      strcpy (REASON, "PSFinitally:ERROR: Could not append PSF.\n");
    }
  else if (NERROR == 332)
    {
      strcpy (REASON, "PSFinitally:ERROR: End-Of-Section mark not found.\n");
    }
  else if (NERROR == 333)
    {
      strcpy (REASON,
	      "IAEAwriteReport:ERROR: unable to write no. of histories to header file.\n");
    }
  else if (NERROR == 334)
    {
      strcpy (REASON,
	      "IAEAwriteReport:ERROR: unable to update header file.\n");
    }
  else if (NERROR == 335)
    {
      strcpy (REASON, "FTLinitally:ERROR: incorrect section header.\n");
    }
  else if (NERROR == 336)
    {
      strcpy (REASON,
	      "FTLinitally:ERROR: Unable to find End-Of-Section mark.\n");
    }
  else if (NERROR == 337)
    {
      strcpy (REASON, "FTLinitally:ERROR: expecting to find ON or OFF.\n");
    }
  else if (NERROR == 338)
    {
      strcpy (REASON, "FTLinitally:ERROR: Too many bins.\n");
    }
  else if (NERROR == 339)
    {
      strcpy (REASON, "FTLinitally:ERROR: End-Of-Section mark not found.\n");
    }
  else if (NERROR == 340)
    {
      strcpy (REASON, "PCSinitally:ERROR: incorrect section header.\n");
    }
  else if (NERROR == 341)
    {
      strcpy (REASON,
	      "PCSinitally:ERROR: Unable to find End-Of-Section mark.\n");
    }
  else if (NERROR == 342)
    {
      strcpy (REASON, "PCSinitally:ERROR: expecting to find ON or OFF.\n");
    }
  else if (NERROR == 343)
    {
      strcpy (REASON, "PCSinitally:ERROR: Too many bins.\n");
    }
  else if (NERROR == 344)
    {
      strcpy (REASON, "PCSinitally:ERROR: End-Of-Section mark not found.\n");
    }
  else if (NERROR == 345)
    {
      strcpy (REASON, "PTSinitally:ERROR: incorrect section header.\n");
    }
  else if (NERROR == 346)
    {
      strcpy (REASON,
	      "PTSinitally:ERROR: Unable to find End-Of-Section mark.\n");
    }
  else if (NERROR == 347)
    {
      strcpy (REASON, "PTSinitally:ERROR: expecting to find ON or OFF.\n");
    }
  else if (NERROR == 348)
    {
      strcpy (REASON, "PTSinitally:ERROR: cannotopen track data file.\n");
    }
  else if (NERROR == 349)
    {
      strcpy (REASON, "PTSinitally:ERROR: End-Of-Section mark not found.\n");
    }
  else if (NERROR == 350)
    {
      strcpy (REASON, "SPDinitally:ERROR: incorrect section header.\n");
    }
  else if (NERROR == 351)
    {
      strcpy (REASON,
	      "SPDinitally:ERROR: Unable to find End-Of-Section mark.\n");
    }
  else if (NERROR == 352)
    {
      strcpy (REASON, "SPDinitally:ERROR: expecting to find ON or OFF.\n");
    }
  else if (NERROR == 353)
    {
      strcpy (REASON, "SPDinitally:ERROR: Invalid entry.\n");
    }
  else if (NERROR == 354)
    {
      strcpy (REASON, "SPDinitally:ERROR: Too many bins.\n");
    }
  else if (NERROR == 355)
    {
      strcpy (REASON, "SPDinitally:ERROR: End-Of-Section mark not found.\n");
    }
  else if (NERROR == 356)
    {
      strcpy (REASON, "PHSinitally:ERROR: incorrect section header.\n");
    }
  else if (NERROR == 357)
    {
      strcpy (REASON,
	      "PHSinitally:ERROR: Unable to find End-Of-Section mark.\n");
    }
  else if (NERROR == 358)
    {
      strcpy (REASON, "PHSinitally:ERROR: expecting to find ON or OFF.\n");
    }
  else if (NERROR == 359)
    {
      strcpy (REASON, "PHSinitally:ERROR: Invalid entry.\n");
    }
  else if (NERROR == 360)
    {
      strcpy (REASON, "PHSinitally:ERROR: Too many bins.\n");
    }
  else if (NERROR == 361)
    {
      strcpy (REASON, "PHSinitally:ERROR: End-Of-Section mark not found.\n");
    }

  else if (NERROR == 362)
    {
      strcpy (REASON, "PIDinitally:ERROR: incorrect section header.\n");
    }
  else if (NERROR == 363)
    {
      strcpy (REASON,
	      "PIDinitally:ERROR: Unable to find End-Of-Section mark.\n");
    }
  else if (NERROR == 364)
    {
      strcpy (REASON, "PIDinitally:ERROR: expecting to find ON or OFF.\n");
    }
  else if (NERROR == 365)
    {
      strcpy (REASON, "PIDinitally:ERROR: Invalid entry.\n");
    }
  else if (NERROR == 366)
    {
      strcpy (REASON,
	      "PIDinitally:ERROR: detection material must be a perfect absorbent; increase absorption energies above einf.\n");
    }
  else if (NERROR == 367)
    {
      strcpy (REASON, "PIDinitally:ERROR: invalid value.\n");
    }
  else if (NERROR == 368)
    {
      strcpy (REASON,
	      "PIDinitally:ERROR: pixel size and no. pixels are both zero.\n");
    }
  else if (NERROR == 369)
    {
      strcpy (REASON, "PIDinitally:ERROR: max no. megapixels exceeded.\n");
    }
  else if (NERROR == 370)
    {
      strcpy (REASON,
	      "PIDinitally:ERROR: emin, emax, nebin invalid values.\n");
    }
  else if (NERROR == 371)
    {
      strcpy (REASON, "PIDinitally:ERROR: invalid mode.\n");
    }
  else if (NERROR == 372)
    {
      strcpy (REASON,
	      "PIDinitally:ERROR: max no. (Mpixels x E bins) exceeded.\n");
    }
  else if (NERROR == 373)
    {
      strcpy (REASON,
	      "PIDinitally:ERROR: c0, c1 parameters cannot be negative.\n");
    }
  else if (NERROR == 374)
    {
      strcpy (REASON,
	      "PIDinitally:ERROR: Matrix format is incompatible with energy discriminating mode.\n");
    }
  else if (NERROR == 375)
    {
      strcpy (REASON, "PIDinitally:ERROR: invalid format.\n");
    }
  else if (NERROR == 376)
    {
      strcpy (REASON, "PIDinitally:ERROR: End-Of-Section mark not found.\n");
    }
  else if (NERROR == 377)
    {
      strcpy (REASON, "PIDinitally:ERROR: not enough memory.\n");
    }
  else if (NERROR == 378)
    {
      strcpy (REASON,
	      "setDetectFrame:ERROR: detection material not found.\n");
    }

  else if (NERROR == 400)
    {
      strcpy (REASON,
	      "PMRDR:ERROR: The input file must begin with the TITLE line.\n");
    }
  else if (NERROR == 401)
    {
      strcpy (REASON, "PMRDR:ERROR: Incorrect particle type.\n");
    }
  else if (NERROR == 402)
    {
      strcpy (REASON,
	      "PMRDR:ERROR: Source energy spectrum. The number of energy bins is too large.\n");
    }
  else if (NERROR == 403)
    {
      strcpy (REASON,
	      "PMRDR:ERROR: The source energy spectrum is not defined.\n");
    }
  else if (NERROR == 404)
    {
      strcpy (REASON, "PMRDR:ERROR: The initial energy E0 is too small.\n");
    }
  else if (NERROR == 405)
    {
      strcpy (REASON, "PMRDR:ERROR: Incorrect body label.\n");
    }
  else if (NERROR == 406)
    {
      strcpy (REASON, "PMRDR:ERROR: THETA must be between 0 and 180 deg.\n");
    }
  else if (NERROR == 407)
    {
      strcpy (REASON, "PMRDR:ERROR: PHI must be between 0 and 360 deg.\n");
    }
  else if (NERROR == 408)
    {
      strcpy (REASON, "PMRDR:ERROR: ALPHA must be between 0 and 180 deg.\n");
    }
  else if (NERROR == 409)
    {
      strcpy (REASON,
	      "PMRDR:ERROR: Inconsistent definition of the primary source.\n");
    }
  else if (NERROR == 410)
    {
      strcpy (REASON, "PMRDR:ERROR: Too many phase-space files.\n");
    }
  else if (NERROR == 411)
    {
      strcpy (REASON, "PMRDR:ERROR: Inconsistent window end points.\n");
    }
  else if (NERROR == 412)
    {
      strcpy (REASON,
	      "PMRDR:ERROR: You have to specify a material file (line MFNAME).\n");
    }
  else if (NERROR == 413)
    {
      strcpy (REASON, "PMRDR:ERROR: Wrong number of materials.\n");
    }
  else if (NERROR == 414)
    {
      strcpy (REASON, "PMRDR:ERROR: Geometry file could not be opened.\n");
    }
  else if (NERROR == 415)
    {
      strcpy (REASON, "PMRDR:ERROR: The PARINP index must be positive.\n");
    }
  else if (NERROR == 416)
    {
      strcpy (REASON, "PMRDR:ERROR: Too many modified parameters.\n");
    }
  else if (NERROR == 417)
    {
      strcpy (REASON, "PMRDR:ERROR: NMATG must be greater than 0.\n");
    }
  else if (NERROR == 418)
    {
      strcpy (REASON, "PMRDR:ERROR: Too many bodies.\n");
    }
  else if (NERROR == 419)
    {
      strcpy (REASON, "PMRDR:ERROR: Too many different materials.\n");
    }
  else if (NERROR == 420)
    {
      strcpy (REASON, "PMRDR:ERROR: Some source bodies are undefined.\n");
    }
  else if (NERROR == 421)
    {
      strcpy (REASON, "PMRDR:ERROR: You have to specify a geometry file.\n");
    }
  else if (NERROR == 422)
    {
      strcpy (REASON, "PMRDR:ERROR: Incorrect body number.\n");
    }
  else if (NERROR == 423)
    {
      strcpy (REASON, "PMRDR:ERROR: Incorrect KB value.\n");
    }
  else if (NERROR == 424)
    {
      strcpy (REASON, "PMRDR:ERROR: Incorrect value of KPAR.\n");
    }
  else if (NERROR == 425)
    {
      strcpy (REASON, "PMRDR:ERROR: Incorrect value of ICOL.\n");
    }
  else if (NERROR == 426)
    {
      strcpy (REASON, "PMRDR:ERROR: Incorrect weight window limits.\n");
    }
  else if (NERROR == 427)
    {
      strcpy (REASON, "PMRDR:ERROR: Incorrect value of IBRSPL.\n");
    }
  else if (NERROR == 428)
    {
      strcpy (REASON,
	      "PMRDR:ERROR: Interaction forcing unactive in this body.\n");
    }
  else if (NERROR == 429)
    {
      strcpy (REASON, "PMRDR:ERROR: Incorrect value of IXRSPL.\n");
    }
  else if (NERROR == 430)
    {
      strcpy (REASON, "PMRDR:ERROR: NBE equal to 0.\n");
    }
  else if (NERROR == 431)
    {
      strcpy (REASON, "PMRDR:ERROR: NBTH equal to 0.\n");
    }
  else if (NERROR == 432)
    {
      strcpy (REASON, "PMRDR:ERROR: Wrong number of PHI bins.\n");
    }
  else if (NERROR == 433)
    {
      strcpy (REASON, "PMRDR:ERROR: Incorrect number of energy bins.\n");
    }
  else if (NERROR == 434)
    {
      strcpy (REASON, "PMRDR:ERROR: Incorrect energy limits.\n");
    }
  else if (NERROR == 435)
    {
      strcpy (REASON, "PMRDR:ERROR: Wrong IPSF value.\n");
    }
  else if (NERROR == 436)
    {
      strcpy (REASON, "PMRDR:ERROR: Wrong IDCUT value.\n");
    }
  else if (NERROR == 437)
    {
      strcpy (REASON,
	      "PMRDR:ERROR: No impact detector has been defined yet.\n");
    }
  else if (NERROR == 438)
    {
      strcpy (REASON,
	      "PMRDR:ERROR: Only one PSF can be generated in a each run.\n");
    }
  else if (NERROR == 439)
    {
      strcpy (REASON, "PMRDR:ERROR: Incorrect number of age bins.\n");
    }
  else if (NERROR == 440)
    {
      strcpy (REASON, "PMRDR:ERROR: Incorrect age limits.\n");
    }
  else if (NERROR == 441)
    {
      strcpy (REASON, "PMRDR:ERROR: Undefined age distribution limits.\n");
    }
  else if (NERROR == 442)
    {
      strcpy (REASON, "PMRDR:ERROR: Incorrect body label.\n");
    }
  else if (NERROR == 443)
    {
      strcpy (REASON,
	      "PMRDR:ERROR: A body cannot be part of two detectors.\n");
    }
  else if (NERROR == 444)
    {
      strcpy (REASON,
	      "PMRDR:ERROR: A void body cannot be part of a detectors.\n");
    }
  else if (NERROR == 445)
    {
      strcpy (REASON, "PMRDR:ERROR: This detector has no active bodies.\n");
    }
  else if (NERROR == 446)
    {
      strcpy (REASON, "PMRDR:ERROR: XU must be greater than XL+1.0E-6.\n");
    }
  else if (NERROR == 447)
    {
      strcpy (REASON, "PMRDR:ERROR: Incorrect keyword.\n");
    }
  else if (NERROR == 448)
    {
      strcpy (REASON, "PMRDR:ERROR: YU must be greater than YL+1.0E-6.\n");
    }
  else if (NERROR == 449)
    {
      strcpy (REASON, "PMRDR:ERROR: ZU must be greater than ZL+1.0E-6.\n");
    }
  else if (NERROR == 450)
    {
      strcpy (REASON, "PMRDR:ERROR: RU must be greater than 1.0E-6.\n");
    }
  else if (NERROR == 451)
    {
      strcpy (REASON,
	      "PMRDR:ERROR: The dump file is corrupted (the TITLE does not match).\n");
    }
  else if (NERROR == 452)
    {
      strcpy (REASON, "PMRDR:ERROR: File could not be opened.\n");
    }
  else if (NERROR == 453)
    {
      strcpy (REASON, "PMRDR:ERROR: The file is empty or corrupted.\n");
    }
  else if (NERROR == 454)
    {
      strcpy (REASON, "ENANG0:ERROR: NBE is too large.\n");
    }
  else if (NERROR == 455)
    {
      strcpy (REASON, "ENANG0:ERROR: NBTH is too large.\n");
    }
  else if (NERROR == 456)
    {
      strcpy (REASON, "ENANG0:ERROR: NBPH is too large.\n");
    }
  else if (NERROR == 457)
    {
      strcpy (REASON,
	      "IMDET0:ERROR: SIMDET: Detector already defined. Detector cannot be defined.\n");
    }
  else if (NERROR == 458)
    {
      strcpy (REASON, "IMDET0:ERROR: SIMDET: Too many detectors.\n");
    }
  else if (NERROR == 459)
    {
      strcpy (REASON, "IMDET0:ERROR: SIMDET: NB is too large.\n");
    }
  else if (NERROR == 460)
    {
      strcpy (REASON,
	      "ENDET0:ERROR: SENDET: Detector already defined. Detector cannot be defined.\n");
    }
  else if (NERROR == 461)
    {
      strcpy (REASON, "ENDET0:ERROR: SENDET: Too many detectors.\n");
    }
  else if (NERROR == 462)
    {
      strcpy (REASON, "ENDET0:ERROR: SENDET: NB is too large.\n");
    }
  else if (NERROR == 463)
    {
      strcpy (REASON, "DOSE0:ERROR: IDOSE should be 1, 2, or 3.\n");
    }
  else if (NERROR == 464)
    {
      strcpy (REASON,
	      "DOSE0:ERROR: SDOSE: NBX must be .GT.0. and .LE.NDXM\n");
    }
  else if (NERROR == 465)
    {
      strcpy (REASON,
	      "DOSE0:ERROR: SDOSE: NBY must be .GT.0. and .LE.NDYM\n");
    }
  else if (NERROR == 466)
    {
      strcpy (REASON,
	      "DOSE0:ERROR: SDOSE: NBZ must be .GT.0. and .LE.NDZM\n");
    }
  else if (NERROR == 1000)
    {
      strcpy (REASON, "");
    }
}
