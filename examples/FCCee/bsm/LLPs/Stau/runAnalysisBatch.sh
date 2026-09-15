#!/bin/bash

theSample=$1

source /afs/desy.de/user/m/mrandria/summer_project/Signal/FCCAnalyses/setup.sh

fccanalysis run analysis_stage1_thomas.py -- --sample $theSample 