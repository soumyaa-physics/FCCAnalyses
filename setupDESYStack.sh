#!/bin/bash

# key4Hep - if you run the script with the "24" argument, you get the old release.
if [[ "$1" == "24" ]]
then 
    source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2024-03-10
else 
    source /cvmfs/sw.hsf.org/key4hep/setup.sh
fi 

# get the directory where this script is located. This should be the FCCAnalysis folder.  
BASE_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )/..

# Add Delphes, assumed to be one folder above FCCAnalysis 
export DELPHES_DIR="${BASE_DIR}/delphes"
export LD_LIBRARY_PATH=${DELPHES_DIR}:$LD_LIBRARY_PATH
export PATH=${DELPHES_DIR}:$PATH 
export CMAKE_PREFIX_PATH=${DELPHES_DIR}:$CMAKE_PREFIX_PATH
export DELPHES_EXTERNALS_TKCOV_INCLUDE_DIR=${DELPHES_DIR}/external/TrackCovariance/

# Add k4SimDelphes, assumed to be one folder above FCCAnalysis
export k4SD_DIR="${BASE_DIR}/k4simdelphes/install/"
export LD_LIBRARY_PATH=${k4SD_DIR}/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${k4SD_DIR}/lib64:$LD_LIBRARY_PATH
export PATH=${k4SD_DIR}/bin:$PATH 
export CMAKE_PREFIX_PATH=${k4SD_DIR}:$CMAKE_PREFIX_PATH

# Add FCCAnalysis, assumed to be the folder we are in now 
export FCCana_DIR="${BASE_DIR}/FCCAnalyses/install/"
export LD_LIBRARY_PATH=${FCCana_DIR}/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${FCCana_DIR}/lib64:$LD_LIBRARY_PATH
export PATH=${FCCana_DIR}/bin:$PATH 
export CMAKE_PREFIX_PATH=${FCCana_DIR}:$CMAKE_PREFIX_PATH
export PYTHONPATH=${FCCana_DIR}/python:$PYTHONPATH

# run this to (re-) compile delphes 
function compileDelphes(){
    cd ${DELPHES_DIR}
    make -j12
    cd - 
}

# run this to configure k4SimDelphes
function confk4SD(){
    mkdir -p ${k4SD_DIR}/../build/
    cd ${k4SD_DIR}/../build/ 
    cmake -DCMAKE_INSTALL_PREFIX=${k4SD_DIR} .. 
    cd - 
}
# run this to (re-)compile k4SimDelphes
function compilek4SD(){
    cmake --build ${k4SD_DIR}/../build -j12 &&  cmake --install ${k4SD_DIR}/../build 
}
# run this to configure fccAnalyses
function confFCC(){
    mkdir -p ${FCCana_DIR}/../build/
    cd ${FCCana_DIR}/../build/ 
    cmake \
        -DWITH_ACTS=OFF \
        -DBUILD_TESTING=OFF \
        -DDELPHES_EXTERNALS_TKCOV_INCLUDE_DIR=$DELPHES_EXTERNALS_TKCOV_INCLUDE_DIR \
        -DCMAKE_INSTALL_PREFIX=${FCCana_DIR} \
        .. 
    cd - 
}
# run this to (re-)compile FCCAnalyses
function compileFCC(){
    cmake --build ${FCCana_DIR}/../build -j12 && cmake --install ${FCCana_DIR}/../build
}