#!/bin/bash

chmod a+x ./src/SimSVGenome.sh
chmod a+x ./src/SVcalling.sh

ln -s $(pwd)/src/SimSVGenome.sh bin/
ln -s $(pwd)/src/SVcalling.sh bin/
ln -s $(pwd)/src/AnnotSV.py bin/
ln -s $(pwd)/src/FilterChimOut.py bin/
ln -s $(pwd)/src/ReadGenome.py bin/
ln -s $(pwd)/src/VerifyFusionGene.py bin/
ln -s $(pwd)/src/VerifySVpred.py bin/

g++ -std=c++11 -o bin/SV2newpos src/SV2newpos.cpp src/GtfTrans.cpp src/SV.cpp src/TRA.cpp src/SimpleSV.cpp -g
g++ -std=c++11 -o bin/GmapSV src/GmapSV.cpp -g
g++ -std=c++11 -o bin/NucmerSV2 src/NucmerSV2.cpp -g
