#!/bin/bash

# key4Hep - if you run the script with the "24" argument, you get the old release.
source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2024-03-10

# get the directory where this script is located. This should be the FCCAnalysis folder. 
BASE_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )/..
export LOCAL_DIR=$(cd $(dirname "${BASH_SOURCE}") && pwd)
# Add Delphes, assumed to be one folder above FCCAnalysis 
export DELPHES_DIR="${BASE_DIR}/delphes"
export LD_LIBRARY_PATH=${DELPHES_DIR}:$LD_LIBRARY_PATH
export PATH=${DELPHES_DIR}:$PATH 
export CMAKE_PREFIX_PATH=${DELPHES_DIR}:$CMAKE_PREFIX_PATH
export DELPHES_EXTERNALS_TKCOV_INCLUDE_DIR=${DELPHES_DIR}/external/TrackCovariance/

# Add FCCAnalysis_pre, assumed to be the folder we are in now 
export FCCana_DIR="${BASE_DIR}/FCCAnalyses_pre/install/"

export LD_LIBRARY_PATH=${FCCana_DIR}/lib:${FCCana_DIR}/lib64:$LD_LIBRARY_PATH
export PATH=${FCCana_DIR}/bin:$PATH
export CMAKE_PREFIX_PATH=${FCCana_DIR}:$CMAKE_PREFIX_PATH
export PYTHONPATH=${FCCana_DIR}/python:$PYTHONPATH
export DELPHES_GEO_PATH=${FCCana_DIR}/../examples/FCCee/bsm/LLPs/Stau


# Make ROOT/Cling use the local FCCAnalyses headers
export ROOT_INCLUDE_PATH="${BASE_DIR}/FCCAnalyses_pre/analyzers/dataframe:${FCCana_DIR}/include:${ROOT_INCLUDE_PATH}"

# run this to (re-) compile delphes 
function compileDelphes(){
    cd ${DELPHES_DIR}
    make -j12
    cd - 
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
    sed -i 's/liblibFCCAnalysis_analysis_example_dictrflx.so/libFCCAnalysis_analysis_example_dictrflx.so/' \
    ${FCCana_DIR}/lib/libFCCAnalysis_analysis_example.rootmap
}