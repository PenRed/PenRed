#!/bin/bash

DIR="$(cd -P "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON=""
PYTHONENV=0

#selection=${1:-"xxx"}
#doseoractivity=${2:-"dose"}
#set1datadir=${3:-"output1"}
#set2datadir=${4:-"output2"}
#labelset1=${5:-"set1"}
#labelset2=${6:-"set2"}
#outputdir=${7:-"Check-${labelset2}-against-${labelset1}"}

function whichpythonisinstalled ()
{
  required="3.6 3.7 3.8 3.9"
  echo ""
  for ver in ${required}
  do
     yesnot=$(command -v python${ver} >/dev/null 2>&1 && echo "1" || echo "0")
     if [[ $yesnot == 1 ]]; then echo "python${ver} is installed"; PYTHON="python${ver}"; else echo "python${ver} is not installed"; fi   
  done
  if [ -z ${PYTHON} ]
  then 
     echo "No required python version is installed."
     echo "Please install one of the required versions: ${required}"
     echo "Aborting ... "
     exit
  else
     echo ""
     echo "${PYTHON} is selected to carry out the checks."
     echo ""
     ispipinstalled=$(${PYTHON} -m pip --version |& grep -q 'No module named pip' && echo '0' || echo '1')
     if [[ $ispipinstalled == 0 ]]
     then
       if [ -d ${DIR}/../env-${PYTHON} ]
       then
          echo -n "${PYTHON} environment found in ${DIR}/../env-${PYTHON}. Activiting ..."
       else
          echo "creating an environment for ${PYTHON}"
          ${PYTHON} -m venv env-${PYTHON}
          echo -n "activiting the environment env-${PYTHON}"
        fi
        source env-${PYTHON}/bin/activate
        echo " ... done"
        echo "If you no longer need it, deactivate it and delete it with rm -rv env-${PYTHON}"
        ${PYTHON} -m pip install --upgrade pip
        PYTHONENV=1
     else
        echo "no environment is needed for ${PYTHON}"
     fi
        echo ""
  fi 
}

function isitinstalled ()
{
  #installed=$(${PYTHON} -m $1 |& grep -q 'is a package' && echo '1' || echo '0')
  aa=($(${PYTHON} -m pip list | grep $1 | awk '{print $1" "$2}'))
  installed=$(if [ -z ${aa[0]} ]; then echo "0"; else if [ "$2" != "*" ]; then if [ "${aa[1]}" == "$2" ]; then echo "1"; else echo "0"; fi; else echo "1"; fi; fi)
  if [[ $installed == 0 ]]
  then 
    #echo "The generation of mhd/raw images requires the python module $1 which is not installed"
    echo "$3"
    case ${silent} in
    "No") 
         echo "Do you wish to install $1 version $2?"
         select yn in "Yes" "No"; do
         case $yn in
           Yes ) if [ "$2" != "*" ]; then ${PYTHON} -m pip install $1==$2; else ${PYTHON} -m pip install $1; fi
                break;;
           No  ) echo "Aborting ..."
                exit
                break;;
         esac
         done
         ;;
    "Yes")
         if [ "$2" != "*" ]; then ${PYTHON} -m pip install $1==$2; else ${PYTHON} -m pip install $1; fi  
         ;;
    esac
  fi
}

function step0 ()
{
  if [ -d "${2}" ]
  then
    echo "Folder ${2} already exists. Deleting ${2}"
    rm -rf ${2}
    echo " "
  else
    echo "Folder ${2} does not exist."
    echo " "
  fi
  #if [ ${1} = "terminate" ];then echo "GoodBye ... "; exit; fi
  if [ ${1} = "terminate" ];then exitscript; exit; fi
}

function step1 ()
{
  echo -n "Step 1: Creating the folder ${outputdir}"
  mkdir -p ${outputdir}/${labelset1}
  mkdir -p ${outputdir}/${labelset2}
  #cd ${outputdir}/PenRed
  echo " ... done"; echo " "
  #if [ ${1} = "terminate" ];then echo "GoodBye ... "; exit; fi
  if [ ${1} = "terminate" ];then exitscript; exit; fi
}

function step2 ()
{
  echo "Step 2: Computing the source activity distribution file SourceActivityDistribution.dat for ${2}"
  if [ ! -d ${outputdir}/${2} ]; then echo " "; echo "Folder ${outputdir} does not exist"; step1 continue; fi
  cd ${outputdir}/${2}
  if [ -f ${DIR}/../${3}/SourceActivityDistribution.dat ]
  then
    echo "Step 2: SourceActivityDistribution found in ${3} for ${2}. We will use it."
    ln -sf ${DIR}/../${3}/SourceActivityDistribution.dat .
  else
    sh ${DIR}/SourceActivityDistribution.sh ${DIR}/../${3}/SpatialSampling.dat ${DIR}/../data/Mhd/${sourcemhd} SourceActivityDistribution.dat "0.0 0.0 0.0"
  fi    
  echo "step 2 done"; echo " "
  #if [ ${1} = "terminate" ];then echo "GoodBye ... "; exit; fi
  if [ ${1} = "terminate" ];then exitscript; exit; fi
  cd ../..
}

function usagestep34 ()
{
  if [ ${1} = "3" ]
  then
    echo "Step 3: flag not recognized."
    echo "Usage: sh Activity-Dose_DistributionChecker.sh step3 [OPTIONS1] [OPTIONS2]"
    echo ""
    echo "Computes the depth dose or depth activity distributions of the patient or the source."
  fi
  if [ ${1} = "4" ]
  then
    echo "Step 4: flag not recognized."
    echo "Usage: sh Activity-Dose_DistributionChecker.sh step4 [OPTIONS1] [OPTIONS2]"
    echo ""
    echo "Generates mhd/raw images of the source activity or patient dose distributions."
  fi
  echo ""
  echo "Options1:"
  echo "terminate  executes the code of the step and exits the script."
  echo "continue   continues the script after executing the step code."
  echo ""
  echo "Options2:"
  if [ ${1} = "3" ]
  then
    echo "activity   computes and compares the ${labelset1} and ${labelset2} depth activity distributions."
    echo "dose       computes and compares the ${labelset1} and ${labelset2} depth dose distributions."
  fi
  if [ ${1} = "4" ]
  then
    echo "activity   computes a mhd/raw image of the source activity distribution."
    echo "dose       computes a mhd/raw image of the patient dose distribution."
  fi
}

function step3 ()
{
  case ${2} in

       activity)
          if [ -f ${DIR}/../${4}/SourceDepthActivity-x.dat ] && [ -f ${DIR}/../${4}/SourceDepthActivity-y.dat ] && [ -f ${DIR}/../${4}/SourceDepthActivity-z.dat ]  
          then
            echo "Step 3: SourceDepthActivity-[xyz].dat files for ${3} found in ${4}. We will use them."
            mkdir -p ${outputdir}/${3}
            cd ${outputdir}/${3}
            ln -sf ${DIR}/../${4}/SourceDepthActivity-*.dat .
          else
            echo "Step 3: Computing depth activity distribution files SourceDepthActivity-[x,y,z].dat for ${3}"
            if [ ! -f ${outputdir}/${3}/SourceActivityDistribution.dat ]; then echo " "; echo "Source activity distribution file does not exist"; step2 continue ${3} ${4}; fi
            cd ${outputdir}/${3}
            sh ${DIR}/ComputeSourceDepthActivity.sh SourceActivityDistribution.dat   ${DIR}/../data/Mhd/${sourcemhd} SourceDepthActivity
          fi
       ;;

       dose)
          if [ -f ${DIR}/../${4}/PatientDepthDose-x.dat ] && [ -f ${DIR}/../${4}/PatientDepthDose-y.dat ] && [ -f ${DIR}/../${4}/PatientDepthDose-z.dat ]
          then
            echo "Step 3: PatientDepthDose-[xyz].dat files for ${3} found in ${4}. We will use them."
            mkdir -p ${outputdir}/${3}
            cd ${outputdir}/${3}
            ln -sf ${DIR}/../${4}/PatientDepthDose-*.dat .
          else
            if [ -z ${PYTHON} ]; then whichpythonisinstalled; fi
            isitinstalled numpy "*" "The calculation of depth dose curves requires the python module numpy which is not installed"
            echo ""
            echo "Step 3: Computing depth dose distribution files PatientDepthDose-[x,y,z].dat for ${3}"
            penreddatadir=${4}
            if [ ! -f ${DIR}/../${penreddatadir}/Dose-th0-spatialDoseDistrib.dat ]; then echo " "; echo "Patient dose distribution file does not exist"; echo "Run the PenRed's pen_main with the input file included in the infiles folder to generate it"; echo "Aborting ..."; exit; fi
            mkdir -p ${outputdir}/${3}
            cd ${outputdir}/${3}
            ${PYTHON} ${DIR}/ComputeDepthDose.py ${DIR}/../${penreddatadir}/Dose-th0-spatialDoseDistrib.dat ${DIR}/../data/Mhd/${patientmhd} PatientDepthDose
          fi  
       ;;

       *)
          usagestep34 3
          exit
       ;;

  esac       
  echo "step 3 done"; echo " "
  #if [ ${1} = "terminate" ];then echo "GoodBye ... "; exit; fi
  if [ ${1} = "terminate" ];then exitscript; exit; fi
  cd ../..
}

function step4 ()
{
  case ${2} in

       activity)
          if [ -f ${DIR}/../${4}/SourceActivity.mhd ] && [ -f ${DIR}/../${4}/SourceActivity.raw ]
          then
            echo "Step 4: SourceActivity.mhd/raw files for ${3} found in ${4}. We will use them."
            mkdir -p ${outputdir}/${3}
            cd ${outputdir}/${3}
            ln -sf ${DIR}/../${4}/SourceActivity.mhd .
            ln -sf ${DIR}/../${4}/SourceActivity.raw .
          else
            if [ -z ${PYTHON} ]; then whichpythonisinstalled; fi
            isitinstalled SimpleITK "*" "The generation of mhd/raw images of the source activity distribution requires the python module SimpleITK which is not installed"
            isitinstalled numpy     "*" "The generation of mhd/raw images of the source activity distribution requires the python module numpy which is not installed"
            echo ""
            echo "Step 4: Generating mhd/raw images of the source activity distribution for ${3}: SourceActivity.mhd/.raw"
            if [ ! -f ${outputdir}/${3}/SourceActivityDistribution.dat ]; then echo " "; echo "Source activity distribution file for ${3} does not exist"; step2 continue ${3} ${4}; fi
            cd ${outputdir}/${3}
            ${PYTHON} ${DIR}/CreateImageFromSourceActivityDistribution.py SourceActivityDistribution.dat "${DIR}/../data/Mhd/${sourcemhd}" SourceActivity 0 
          fi  
       ;;

       dose)
          if [ -f ${DIR}/../${4}/PatientDose.mhd ] && [ -f ${DIR}/../${4}/PatientDose.raw ]
          then
            echo "Step 4: PatientDose.mhd/raw files for ${3} found in ${4}. We will use them."
            mkdir -p ${outputdir}/${3}
            cd ${outputdir}/${3}
            ln -sf ${DIR}/../${4}/PatientDose.mhd .
            ln -sf ${DIR}/../${4}/PatientDose.raw .
          else
            if [ -z ${PYTHON} ]; then whichpythonisinstalled; fi
            isitinstalled SimpleITK "*" "The generation of mhd/raw images of the patient dose distribution requires the python module SimpleITK which is not installed"
            isitinstalled numpy     "*" "The generation of mhd/raw images of the patient dose distribution requires the python module numpy which is not installed"
            echo ""
            echo "Step 4: Generating mhd/raw images of the patient dose distribution for ${3}: PatientDose.mhd/.raw"
            penreddatadir=${4}
            if [ ! -f ${DIR}/../${penreddatadir}/Dose-th0-spatialDoseDistrib.dat ]; then echo " "; echo "Patient dose distribution file for ${3} does not exist"; echo "Run PenRed's pen_main with the input file included in the infiles folder to generate it"; echo "Aborting ..."; exit; fi
            mkdir -p ${outputdir}/${3}
            cd ${outputdir}/${3}
            ${PYTHON} ${DIR}/CreateImageFromSpatialDoseDistribution.py ${DIR}/../${penreddatadir}/Dose-th0-spatialDoseDistrib.dat   "${DIR}/../data/Mhd/${patientmhd}" PatientDose 0
          fi  
       ;;

       *)
          usagestep34 4
          exit
       ;;

  esac    
  echo "step 4 done"; echo " "
  #if [ ${1} = "terminate" ];then echo "GoodBye ... "; exit; fi
  if [ ${1} = "terminate" ];then exitscript; exit; fi
  cd ../..
}

function check1 ()
{
  case ${2} in

       activity)
          echo "Check 1: comparision of the ${labelset1} and ${labelset2} depth activity distributions of the source"
          if [ ! -f ${outputdir}/${labelset1}/SourceDepthActivity-x.dat ]; then echo " "; echo "Source depth activity distribution files for ${labelset1} do not exist"; step3 continue activity ${labelset1} ${set1datadir}; fi
          if [ ! -f ${outputdir}/${labelset2}/SourceDepthActivity-x.dat ]; then echo " "; echo "Source depth activity distribution files for ${labelset2} do not exist"; step3 continue activity ${labelset2} ${set2datadir}; fi
          mkdir -p ${outputdir}/DepthActivityDistribution; cd ${outputdir}/DepthActivityDistribution
          ln -sf ../${labelset1} .
          ln -sf ../${labelset2} .
          sh ${DIR}/Gnuplot/plotDepthActivity.sh ${DIR}/../data/Mhd/${sourcemhd} SourceDepthActivity SourceDepthActivity ${labelset1} ${labelset2} ${labelset1} ${labelset2}
          echo "Figure png files for comparison of the depth activity created."
       ;;

       dose)
          echo "Check 1: comparision of the ${labelset1} and ${labelset2} depth dose distributions of the patient"
          if [ ! -f ${outputdir}/${labelset1}/PatientDepthDose-x.dat ]; then echo " "; echo "Patient depth dose distribution files for ${labelset1} do not exist"; step3 continue dose ${labelset1} ${set1datadir}; fi
          if [ ! -f ${outputdir}/${labelset2}/PatientDepthDose-x.dat ]; then echo " "; echo "Patient depth dose distribution files for ${labelset2} do not exist"; step3 continue dose ${labelset2} ${set2datadir}; fi
          mkdir -p ${outputdir}/DepthDoseDistribution; cd ${outputdir}/DepthDoseDistribution
          ln -sf ../${labelset1} .
          ln -sf ../${labelset2} .
          sh ${DIR}/Gnuplot/plotDepthDose.sh ${DIR}/../data/Mhd/${patientmhd} PatientDepthDose PatientDepthDose ${labelset1} ${labelset2} ${labelset1} ${labelset2}
          echo "Figure png files for comparison of the depth dose created."
       ;;

       *)
          usage
          exit
       ;;       

  esac
  case ${silent} in
  "No") 
       echo "Do you wish to see these curves with okular?"
       select yn in "Yes" "No"; do
         case $yn in
           Yes ) for ifile in `ls -1 *.png`
                 do
                    okular ${ifile}
                 done; 
                 break;;
           No ) break;;
         esac
       done
       ;;
  "Yes") ;;
  esac
  echo "check 1 done"; echo " " 
  #if [ ${1} = "terminate" ];then echo "GoodBye ... "; exit; fi
  if [ ${1} = "terminate" ];then exitscript; exit; fi
  cd ../..
}

function check2 ()
{

  case ${2} in

       activity)
          echo "Check 2: computation of the gamma index between the Gate and PenRed spatial activity distributions"
          mhdfile="SourceActivity.mhd"
          if [ ! -f ${outputdir}/${labelset1}/${mhdfile} ]; then echo " "; echo "Mhd/raw image of the source activity distribution for ${labelset1} does not exist"; step4 continue activity ${labelset1} ${set1datadir}; fi
          if [ ! -f ${outputdir}/${labelset2}/${mhdfile} ]; then echo " "; echo "Mhd/raw image of the source activity distribution for ${labelset2} does not exist"; step4 continue activity ${labelset2} ${set2datadir}; fi
          mkdir -p ${outputdir}/GammaIndex-SourceActivity; cd ${outputdir}/GammaIndex-SourceActivity
          dd[0]=1; dd[1]=0.5
          ddunit[0]="%"; ddunit[1]="%"
          dta[0]=1.0; dta[1]=0.5
          T[0]=1.0; T[1]=0.5
          logfile="SourceActivity-gamma.log"
       ;;

       dose)
          echo "Check 2: computation of the gamma index between the Gate and PenRed spatial dose distributions"
          mhdfile="PatientDose.mhd"
          if [ ! -f ${outputdir}/${labelset1}/${mhdfile} ]; then echo " "; echo "Mhd/raw image of the patient dose distribution for ${labelset1} does not exist"; step4 continue dose ${labelset1} ${set1datadir}; fi
          if [ ! -f ${outputdir}/${labelset2}/${mhdfile} ]; then echo " "; echo "Mhd/raw image of the patient dose distribution for ${labelset2} does not exist"; step4 continue dose ${labelset2} ${set2datadir}; fi
          mkdir -p ${outputdir}/GammaIndex-PatientDose; cd ${outputdir}/GammaIndex-PatientDose
          dd[0]=1; dd[1]=0.5; dd[2]=0.5
          ddunit[0]="%"; ddunit[1]="%"; ddunit[2]="%"
          dta[0]=1.0; dta[1]=0.5; dta[2]=0.5
          T[0]=1.0; T[1]=0.5; T[2]=0.0001
          logfile="PatientDose-gamma.log"
       ;;        
    
       *)
          usage
          exit
       ;;       

  esac
  if [ -z ${PYTHON} ]; then whichpythonisinstalled; fi
  isitinstalled tqdm "4.50.2" "The computation of the gamma index both of dose and activity distributions requires the python module tqdm 4.50.2 which is not installed"
  isitinstalled psutil    "*" "The computation of the gamma index both of dose and activity distributions requires the python module psutil which is not installed"
  isitinstalled gatetools "*" "The computation of the gamma index both of dose and activity distributions requires the python module gatetools which is not installed"
  isitinstalled SimpleITK "*" "The computation of the gamma index both of dose and activity distributions requires the python module SimpleITK which is not installed"
  isitinstalled numpy     "*" "The computation of the gamma index both of dose and activity distributions requires the python module numpy which is not installed"
  echo ""
  for ipar in `seq 0 $((${#dd[@]}-1))`
  do
    echo "computing gamma index between reference ${labelset1} and target ${labelset2} using dd=${dd[$ipar]} ddunit=${ddunit[$ipar]} dta=${dta[$ipar]} T=${T[$ipar]}"
    gammaindexparam="${dd[$ipar]}-${dta[$ipar]}-${T[$ipar]}"
    if [ ${silent} = "No" ]
    then
       gt_gamma_index ../${labelset1}/${mhdfile} ../${labelset2}/${mhdfile} -q -o aux-gamma.mhd --dd ${dd[$ipar]} --dta ${dta[$ipar]} -u "${ddunit[$ipar]}" -T ${T[$ipar]} -v --logfile ${logfile}
    else   
       gt_gamma_index ../${labelset1}/${mhdfile} ../${labelset2}/${mhdfile} -q -o aux-gamma.mhd --dd ${dd[$ipar]} --dta ${dta[$ipar]} -u "${ddunit[$ipar]}" -T ${T[$ipar]} --logfile ${logfile}
    fi   
    resultsfile="GammaIndex-${gammaindexparam}.res"
    ${PYTHON} ${DIR}/CalculateGammaIndex.py aux-gamma.mhd ${dd[$ipar]} ${ddunit[$ipar]} ${dta[$ipar]} ${T[$ipar]} &> ${resultsfile}
    more ${resultsfile}
    echo ""
    rm -f aux-gamma.mhd; rm -f aux-gamma.raw
    mkdir -p GammaIndex-${gammaindexparam}; mv ${resultsfile} ${logfile} GammaIndex-${gammaindexparam}
  done
  echo "check 2 done"; echo " "
  #if [ ${1} = "terminate" ];then echo "GoodBye ... "; exit; fi
  if [ ${1} = "terminate" ];then exitscript; exit; fi
  cd ../..
}

function zsliceSelection ()
{
  echo "" | awk -v nZ="${1}" -v sZ="${2}" 'BEGIN{sl=6; ns=4; ps1=0.15; ps2=0.05}{lower=int(nZ*ps1); upper=nZ-int(nZ*ps2); dd=int((upper-lower)/(sl*ns-2)); ic=-1; string=""; for(ii=lower; ii<=upper; ii+=dd){ic++; if((ic-ic%sl)/sl==sZ){str=sprintf("%d",ii); string=string" "str;}}}END{print string}'
}

function check3 ()
{

  case ${2} in

       activity)
          if [ -z ${PYTHON} ]; then whichpythonisinstalled; fi
          isitinstalled SimpleITK "*" "The generation of colormap 2D images of the source activty distribution requires the python module SimpleITK which is not installed"
          isitinstalled numpy     "*" "The generation of colormap 2D images of the source activty distribution requires the python module numpy which is not installed"
          echo ""
          echo "Check 3: colormap 2D images of the source activity distribution in selected Z slices"
          if [ ! -f ${outputdir}/${labelset1}/SourceActivity.mhd ]; then echo " "; echo "Mhd/raw image of the source activity distribution for ${labelset1} does not exist"; step4 continue activity ${labelset1} ${set1datadir}; fi
          if [ ! -f ${outputdir}/${labelset2}/SourceActivity.mhd ]; then echo " "; echo "Mhd/raw image of the source activity distribution for ${labelset2} does not exist"; step4 continue activity ${labelset2} ${set2datadir}; fi
          mkdir -p ${outputdir}/Color2DMaps-SourceActivity; cd ${outputdir}/Color2DMaps-SourceActivity
          echo -n "merging ${labelset1} and ${labelset2} source activity distributions"
          ${PYTHON} ${DIR}/MergeSourceActivityDistributions-Images.py ../${labelset1}/SourceActivity.mhd ../${labelset2}/SourceActivity.mhd merged.dat
          echo " ... done"
          echo -n "creating both png single colormap plots and pdf comparison files"
          mkdir -p SingleColorMaps/${labelset1}; mkdir -p SingleColorMaps/${labelset2}; mkdir -p ComparisonColorMaps
          datamhdfile="${sourcemhd}"
          titleSet1="Source activity distr. from ${labelset1}"
          titleSet2="Source activity distr. from ${labelset2}"
          declare -a zslicelist
          nZ=$(more ../${labelset1}/SourceActivity.mhd | awk '{if(/^DimSize/){print $5}}')

          zslicelist[0]="$(zsliceSelection $nZ 0)"
          zslicelist[1]="$(zsliceSelection $nZ 1)"
          zslicelist[2]="$(zsliceSelection $nZ 2)"
          zslicelist[3]="$(zsliceSelection $nZ 3)"

          #zslicelist[0]="30 35 40 45 50 55"
          #zslicelist[1]="60 65 70 75 80 85" 
          #zslicelist[2]="90 95 101 120 140 160"
          #zslicelist[3]="165 170 175 180 185 190"
       ;;

       dose)
          if [ -z ${PYTHON} ]; then whichpythonisinstalled; fi
          isitinstalled SimpleITK "*" "The generation of colormap 2D images of the patient dose distribution requires the python module SimpleITK which is not installed"
          isitinstalled numpy     "*" "The generation of colormap 2D images of the patient dose distribution requires the python module numpy which is not installed"
          echo ""
          echo "Check 3: colormap 2D images of the patient dose distribution in selected Z slices"
          if [ ! -f ${outputdir}/${labelset1}/PatientDose.mhd ]; then echo " "; echo "Mhd/raw image of the patient dose distribution for ${labelset1} does not exist"; step4 continue dose ${labelset1} ${set1datadir}; fi
          if [ ! -f ${outputdir}/${labelset2}/PatientDose.mhd ]; then echo " "; echo "Mhd/raw image of the patient dose distribution for ${labelset2} does not exist"; step4 continue dose ${labelset2} ${set2datadir}; fi
          mkdir -p ${outputdir}/Color2DMaps-PatientDose; cd ${outputdir}/Color2DMaps-PatientDose
          echo -n "merging ${labelset1} and ${labelset2} patient dose spatial distributions"
          if [ ! -f merged.dat ]
          then
            ${PYTHON} ${DIR}/MergeSourceActivityDistributions-Images.py ../${labelset1}/PatientDose.mhd ../${labelset2}/PatientDose.mhd merged.dat
            echo " ... done"
          else
            echo ""
            echo "merged patient dose distribution file found in ${outputdir}/Color2DMaps-PatientDose. We will use it."
          fi  
          echo -n "creating both png single colormap plots and pdf comparison files"
          mkdir -p SingleColorMaps/${labelset1}; mkdir -p SingleColorMaps/${labelset2}; mkdir -p ComparisonColorMaps
          datamhdfile="${patientmhd}"
          titleSet1="Patient dose distr. from ${labelset1}"
          titleSet2="Patient dose distr. from ${labelset2}"
          declare -a zslicelist
          nZ=$(more ../${labelset1}/PatientDose.mhd | awk '{if(/^DimSize/){print $5}}')

          zslicelist[0]="$(zsliceSelection $nZ 0)"
          zslicelist[1]="$(zsliceSelection $nZ 1)"
          zslicelist[2]="$(zsliceSelection $nZ 2)"
          zslicelist[3]="$(zsliceSelection $nZ 3)"

          #zslicelist[0]="18 22 32 42 47 50"
          #zslicelist[1]="55 60 65 70 75 80" 
          #zslicelist[2]="90 100 110 120 140 160"
          #zslicelist[3]="165 170 175 180 185 190"
       ;;

       *)
          usage
          exit
       ;;       

  esac
  for izsl in "${zslicelist[@]}"
  do
     zslicestring=`echo ${izsl} | sed 's/\ /-/g'`
     sh ${DIR}/Gnuplot/plotImage.sh ${DIR}/../data/Mhd/${datamhdfile} "${izsl}" "${titleSet1}" "${titleSet2}" "${labelset1}" "${labelset2}" &> /dev/null
     mv image.pdf ComparisonColorMaps/image-${zslicestring}.pdf   
  done 
  rm -f merged.dat
  mv ${labelset1}-*.png SingleColorMaps/${labelset1}
  mv ${labelset2}-*.png SingleColorMaps/${labelset2}
  echo " ... done"
  case ${silent} in
  "No") 
       echo "Do you wish to see the comparison pdf files with okular?"
       select yn in "Yes" "No"; do
         case $yn in
           Yes ) for ifile in `ls -tr ComparisonColorMaps/*.pdf`
                 do
                    okular ${ifile}
                 done; 
                 break;;
           No ) break;;
         esac
       done
       ;;
  "Yes") ;;
  esac
  echo "check 3 done"; echo " " 
  #if [ ${1} = "terminate" ];then echo "GoodBye ... "; exit; fi
  if [ ${1} = "terminate" ];then exitscript; exit; fi
  cd ../..
}

function exitscript ()
{

  if [ -d ${DIR}/../env-${PYTHON} ]
  then
    echo "${PYTHON} environment found in ${DIR}/../env-${PYTHON}."
    echo "If you no longer need it, deactivate it and delete it with rm -rv env-${PYTHON}"
    echo ""
  fi

  if [[ $PYTHONENV == 1 ]]
  then
        echo -n "deactiviting the env"
        deactivate
        echo " ... done"
        PYTHONENV=0
        echo ""
  fi

  echo "Goodbye ...."

}

function usageOLD ()
{
  echo "Usage: sh Activity-Dose_DistributionChecker.sh [OPTIONS] [QUANTITY] [REFERENCEDIR] [TARGETDIR] [OUTPUTDIR]"
  echo ""
  echo "Checks ${labelset2} results against ${labelset1}"
  echo ""
  echo "OPTIONS:"
  echo "check1     Compares ${labelset2} source depth activity distribution or patient depth dose distribution to the ${labelset1} ones."
  echo "check2     Computes the gamma index between the ${labelset1} and ${labelset2} source depth activity or patient depth dose distributions."
  echo "check3     Computes and compares colormap 2D images of ${labelset1} and ${labelset2} source activity or patient dose distributions in selected Z slices"
  echo "allchecks  Performs all checks."
  echo "help       Show this message and exit."
  echo ""
  echo "QUANTITY:"
  echo "activity   Compares source activity distributions."
  echo "dose       Compares patient dose distributions."
  echo ""
  echo "REFERENCEDIR: Folder containing the ${labelset1} reference simulation data."
  echo "              Optional. Default: output1"
  echo "TARGETDIR:    Folder containing the ${labelset2} target simulation data."
  echo "              Optional. Default: output2"
  echo "OUTPUTDIR:    Folder in which the analisis files will be saved."
  echo "              Optional. Default: Check-${labelset2}-against-${labelset1}"
  echo ""
}

function usage ()
{
  echo "Usage: sh Activity-Dose_DistributionChecker.sh -c [CHECKS] -q [QUANTITY] -r [REFERENCEDIR] -t [TARGETDIR] -R [LABELREFERENCE] -T [LABELTARGET] -o [OUTPUTDIR] -mM [METADATA] -s [SILENT]"
  echo ""
  echo "Checks ${labelset2} simulation results against ${labelset1} ones"
  echo ""
  echo "CHECKS:         default = none"
  echo "check1          Compares ${labelset2} source depth activity distribution or patient depth dose distribution to the ${labelset1} ones."
  echo "check2          Computes the gamma index between the ${labelset1} and ${labelset2} source depth activity or patient depth dose distributions."
  echo "check3          Computes and compares colormap 2D images of ${labelset1} and ${labelset2} source activity or patient dose distributions in selected Z slices"
  echo "allchecks|all   Performs all checks."
  echo "help            Show this message and exit."
  echo ""
  echo "QUANTITY:       default = dose"
  echo "activity        Compares source activity distributions."
  echo "dose            Compares patient dose distributions. Default."
  echo "all             Compares both activity and dose distributions."
  echo ""
  echo "REFERENCEDIR:   default = output1"   
  echo "                Folder containing the ${labelset1} reference simulation data."
  echo ""
  echo "LABELREFERENCE: default = set1"
  echo "                Label/identification for the reference simulation data."
  echo ""
  echo "TARGETDIR:      default = output2"
  echo "                Folder containing the ${labelset2} target simulation data."
  echo ""
  echo "LABELTARGET:    default = set2"
  echo "                Label/identification for the target simulation data."
  echo ""
  echo "OUTPUTDIR:      default = Check-LABELTARGET-against-LABELREFERENCE"
  echo "                Folder in which the analisis files will be saved."
  echo ""
  echo "METADATA:       default = for -m, Patient.mhd; for -M, Source.mhd"
  echo "                Files containing the metadata of the patient and source images: number of dimensions, offset,"
  echo "                spacings and dimension sizes among others. These files must be in the folder data/Mhd" 
  echo ""
  echo "SILENT:         default = verbose"
  echo "                Sets verbose o silent screen output"
  echo ""
}

selection="xxx"
doseoractivity="dose"
set1datadir="output1"
set2datadir="output2"
labelset1="set1"
labelset2="set2"
outputdir="Check-${labelset2}-against-${labelset1}"
otriggered=0
patientmhd="Patient.mhd"
sourcemhd="Source.mhd"
silent="No"

while getopts "c:q:r:t:R:T:o:m:M:s" option
do
   case "${option}" in
     c) selection="${OPTARG}";;
     q) doseoractivity="${OPTARG}";;
     r) set1datadir="${OPTARG}";;
     t) set2datadir="${OPTARG}";;
     R) labelset1="${OPTARG}";;
     T) labelset2="${OPTARG}";;
     o) outputdir="${OPTARG}"; otriggered=1;;
     m) patientmhd="${OPTARG}";;
     M) sourcemhd="${OPTARG}";;
     s) silent="Yes";;
    \?) usage; exit;;
   esac
done

if [[ ${otriggered} = 0 ]]; then outputdir="Check-${labelset2}-against-${labelset1}"; fi

#echo "Your selection: -c ${selection} -q ${doseoractivity} -r ${set1datadir} -t ${set2datadir} -R ${labelset1} -T ${labelset2} -o ${outputdir} -s ${silent}"

if [ ${silent} = "Yes" ]; then exec 3>&1 &>/dev/null; fi
#command_you_want_to_see >&3
  
case $selection in

     help)
        usage
     ;;
     delete)
        step0 terminate ${outputdir}
     ;;
     check1) 
        check1 terminate ${doseoractivity}
     ;;
     check2) 
        check2 terminate ${doseoractivity}
     ;;
     check3) 
        check3 terminate ${doseoractivity}
     ;;
     allchecks|all)
        case ${doseoractivity} in
          activity|dose)
             echo "Executing all checks for ${doseoractivity} distribution"
             check1 continue ${doseoractivity}
             check2 continue ${doseoractivity}
             check3 continue ${doseoractivity}
          ;;
          all)
             echo "Executing all checks for both activity and dose distributions"
             check1 continue activity
             check1 continue dose
             check2 continue activity
             check2 continue dose
             check3 continue activity
             check3 continue dose
          ;;
          *)
             usage
          ;;
        esac
     ;;       
     *)
        usage
     ;;

esac

if [ -d ${DIR}/../env-${PYTHON} ]
then
    echo "${PYTHON} environment found in ${DIR}/../env-${PYTHON}."
    echo "If you no longer need it, deactivate it and delete it with rm -rv env-${PYTHON}"
    echo ""
fi

if [[ $PYTHONENV == 1 ]]
then
        echo -n "deactiviting the env"
        deactivate
        echo " ... done"
        PYTHONENV=0
        echo ""
fi

echo "Goodbye ...."
