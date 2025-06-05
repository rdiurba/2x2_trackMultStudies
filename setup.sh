source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh

setup duneanaobj v03_06_01b -q e20:prof
setup dk2nugenie   v01_10_01k -q debug:e20
setup cmake v3_27_4
# edep-sim needs to have the GEANT bin directory in the path
export PATH=$PATH:$GEANT4_FQ_DIR/bin

# shut up ROOT include errors
export ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH:$GENIE_INC/GENIE

# nusystematics paths
#export NUSYST=${PWD}/nusystematics
#export LD_LIBRARY_PATH=${NUSYST}/build/Linux/lib:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=${NUSYST}/build/nusystematics/artless:$LD_LIBRARY_PATH
#export FHICL_FILE_PATH=${NUSYST}/nusystematics/fcl:$FHICL_FILE_PATH

# Add pyGeoEff to pythonpath
export PYTHONPATH=${PYTHONPATH}:${PWD}/DUNE_ND_GeoEff/lib/

# duneananobj needs to be in the libs too
export LD_LIBRARY_PATH=${DUNEANAOBJ_LIB}:$LD_LIBRARY_PATH

# finally, add our lib & bin to the paths
mydir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
export LD_LIBRARY_PATH=$mydir/lib:$LD_LIBRARY_PATH
export PATH=$mydir/bin:$PATH


export ROOUNFOLD=/global/homes/r/rdiurba/cafAna/RooUnfold
export LD_LIBRARY_PATH=$ROOUNFOLD:$LD_LIBRARY_PATH
export C_INCLUDE_PATH="$ROOUNFOLD:$C_INCLUDEPATH"
export CPLUS_INCLUDE_PATH=$ROOUNFOLD:$CPLUS_INCLUDE_PATH
export PATH="$ROOUNFOLD:$PATH"

export LD_LIBRARY_PATH=$ROOUNFOLD/src:$LD_LIBRARY_PATH
export C_INCLUDE_PATH="$ROOUNFOLD/src:$C_INCLUDEPATH"
export CPLUS_INCLUDE_PATH=$ROOUNFOLD/src:$CPLUS_INCLUDE_PATH
export PATH="$ROOUNFOLD/src:$PATH"

