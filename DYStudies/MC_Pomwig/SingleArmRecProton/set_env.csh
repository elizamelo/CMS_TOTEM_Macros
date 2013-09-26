stty erase "^?"
#function setenv() { export "$1=$2"; }
setenv LC_ALL en_US
#cd /afs/cern.ch/exp/totem/scratch/hubert/cmssw_code/CMSSW_4_2_4_SVN_4/CMSSW_4_2_4/src
cd /afs/cern.ch/exp/totem/scratch/hubert/cmssw_code/CMSSW_4_2_4_SVN_6/CMSSW_4_2_4/src/
eval `scram runtime -csh`
cd -
