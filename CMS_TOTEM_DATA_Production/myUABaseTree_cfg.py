from UATree.UABaseTree.UABaseTree_forward_cfg import process
process.source.fileNames = ['file:/storage2/sfonseca/CMSSW/TOTEM/CMSSW_5_3_2_patch4/src/test/test.root']
process.maxEvents.input = 1000
process.GlobalTag.globaltag = 'FT_R_53_V10::All'

process.uabasetree.outputfilename = "UABaseTree_CMS-TOTEM.root"
