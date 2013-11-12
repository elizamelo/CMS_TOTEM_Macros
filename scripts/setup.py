import sys
import optparse
from subprocess import Popen,PIPE
#import TOTEM_CMS_Combine_Files
from TOTEM_CMS_Combine_Files import*

## Parameters
castor_dir = "/castor/cern.ch/totem/offline/CMSTOTEM/TotemNtuples/HighBeta/8369"
local_dir = "/afs/cern.ch/work/s/sfonseca/CMS_TOTEM_Merge/CMS_LP_ExclEGMU_Run2012C-22Jan2013-v1_RECO_Run_198902"
type ="root"
##################################################################
def combine(castor_dir,local_dir,type):
     from subprocess import call
     files_ = TOTEM_CMS_Combine_Files(castor_dir,local_dir,type)
     return files_
##########################################################################
#directory where totem ntuples are stored, please put "rfio:directoryName" if it is on CASTOR
totemDirectory = "rfio:" + castor_dir 

#directory where cms ntuples are stored, "rfio:directoryName" if it's on CASTOR
cmsDirectory = "file:" + local_dir

#the map contains the runs that are going to be merged
#key - merged file name
#value - [totemNtuple, cmsNtuple]

files = combine(castor_dir,local_dir,type)

#if the output directory is on CASTOR please set outputCastor = True
outputCastor = True
#a place where merged files will be stored
outputDirectory = "/castor/cern.ch/user/s/sfonseca/8369_TOTEM_ALL_MERGED_13JUN"

#the path to the compiled program
compiledProgramPath = "/afs/cern.ch/work/s/sfonseca/CMS_TOTEM_Merge/CMSTotem/Merge/mergeNtuples"


#the queue which will be used on lxbatch
queueName = "1nd"

print "========================="
print "----> Combine completed."
print "========================="



