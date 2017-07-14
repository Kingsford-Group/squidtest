#!/bin/bash

cp scripts/SimSVGenome.sh bin/
cp scripts/SVcalling.sh bin/
cp scripts/extractSTAR.sh bin/
chmod a+x ./bin/SimSVGenome.sh
chmod a+x ./bin/SVcalling.sh
chmod a+x ./bin/extractSTAR.sh
g++ -std=c++11 -o bin/SV2newpos scripts/SV2newpos.cpp scripts/GtfTrans.cpp scripts/SV.cpp scripts/TRA.cpp scripts/SimpleSV.cpp -g
g++ -std=c++11 -o bin/GmapSV scripts/GmapSV.cpp -g
g++ -std=c++11 -o bin/NucmerSV2 scripts/NucmerSV2.cpp -g
g++ -std=c++11 -o bin/mergesplit scripts/mergesplit.cpp -I ${BamtoolsDir}/include/ -L ${BamtoolsDir}/lib/ -lbamtools -lz -g -Wl,-rpath,${BamtoolsDir}/lib
