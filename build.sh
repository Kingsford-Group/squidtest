#!/bin/bash

args=("$@")
CPPBOOST=""
if [[ ${args[0]} == "--BOOST" ]] && ((${#args}>1)); then
	CPPBOOST="-I"${args[1]}
fi

mkdir -p bin

chmod a+x ./src/SimSVGenome.sh
chmod a+x ./src/SVcalling.sh

if [[ ! -L $(pwd)/bin/SimSVGenome.sh ]] || [[ ! -e $(pwd)/bin/SimSVGenome.sh ]]; then
	ln -sf $(pwd)/src/SimSVGenome.sh bin/
fi
if [[ ! -L $(pwd)/bin/SVcalling.sh ]] || [[ ! -e $(pwd)/bin/SVcalling.sh ]]; then
	ln -sf $(pwd)/src/SVcalling.sh bin/
fi
if [[ ! -L $(pwd)/bin/AnnotSV.py ]] || [[ ! -e $(pwd)/bin/AnnotSV.py ]]; then
	ln -sf $(pwd)/src/AnnotSV.py bin/
fi
if [[ ! -L $(pwd)/bin/FilterChimOut.py ]] || [[ ! -e $(pwd)/bin/FilterChimOut.py ]]; then
	ln -sf $(pwd)/src/FilterChimOut.py bin/
fi
if [[ ! -L $(pwd)/bin/ReadGenome.py ]] || [[ ! -e $(pwd)/bin/ReadGenome.py ]]; then
	ln -sf $(pwd)/src/ReadGenome.py bin/
fi
if [[ ! -L $(pwd)/bin/VerifyFusionGene.py ]] || [[ ! -e $(pwd)/bin/VerifyFusionGene.py ]]; then
	ln -sf $(pwd)/src/VerifyFusionGene.py bin/
fi
if [[ ! -L $(pwd)/bin/VerifySVpred.py ]] || [[ ! -e $(pwd)/bin/VerifySVpred.py ]]; then
	ln -sf $(pwd)/src/VerifySVpred.py bin/
fi

g++ -std=c++11 -o bin/SV2newpos src/SV2newpos.cpp src/GtfTrans.cpp src/SV.cpp src/TRA.cpp src/SimpleSV.cpp ${CPPBOOST} -g
g++ -std=c++11 -o bin/GmapSV src/GmapSV.cpp ${CPPBOOST} -g
g++ -std=c++11 -o bin/NucmerSV2 src/NucmerSV2.cpp ${CPPBOOST} -g
