#!/bin/bash

#This script downloads all the available nucleide data in PENNUC format

#Get nucleide names
curl -s http://www.nucleide.org/DDEP_WG/DDEPdata.htm | grep '<td><font size="-1"><b>' | grep '</b></font></td>' | cut -d ">" -f 4 | cut -d "<" -f 1 &> nucleideNames.txt

#Download files
rm -r data
mkdir data
cat nucleideNames.txt | while read line; do wget -nv -O data/${line} http://www.nucleide.org/DDEP_WG/Nuclides/${line}.PenNuc.txt; done
