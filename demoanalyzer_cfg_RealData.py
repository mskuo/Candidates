import FWCore.ParameterSet.Config as cms
from RecoMuon.TrackingTools.MuonServiceProxy_cff import *
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes
process = cms.Process("Demo")

process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

#globaltag for 2012 collision data
process.GlobalTag.globaltag = 'FT53_V21A_AN6::All'

# intialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('Demo')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
        limit = cms.untracked.int32(-1)
        )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
# **********************************************************************
# set the maximum number of events to be processed                     *
#    this number (argument of int32) is to be modified by the user     *
#    according to need and wish                                        *
#    default is preset to -1 (all events)                              *
# **********************************************************************
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000000))


# set the number of events to be skipped (if any) at end of file below

# define JSON file for 2012 data
goodJSON = '/home/menghsiu/CMSSW_5_3_32/src/Demo/DemoAnalyzer/JSON/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt'

myLumis = LumiList.LumiList(filename = goodJSON).getCMSSWString().split(',')

# ****************************************************************************
# define the input data set here by inserting the appropriate .txt file list *
# ****************************************************************************
import FWCore.Utilities.FileUtils as FileUtils
#
# ********************************************************************
# load the data set                                                  * 
# this example uses one subset of the relevant 2012 DoubleMu dataset *
# ********************************************************************
#
# use the following if you want to run over a full index file
files2012data = FileUtils.loadListFromFile ('/home/menghsiu/CMSSW_5_3_32/src/Demo/DemoAnalyzer/datasets/CMS_Run2012B_DoubleMuParked_AOD_22Jan2013-v1_20000_file_index.txt')
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(*files2012data    
    )    
)

# use the file in your local
#process.source = cms.Source("PoolSource",
#	fileNames = cms.untracked.vstring ( 'file:/home/menghsiu/CMSSW_5_3_32/src/Demo/DemoAnalyzer/001E739E-9268-E211-A1E1-00259073E32A.root' ) 
#	)

# apply JSON file
#   (needs to be placed *after* the process.source input file definition!)
process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
process.source.lumisToProcess.extend(myLumis)

# *************************************************
# number of events to be skipped (0 by default)   *
# *************************************************
process.source.skipEvents = cms.untracked.uint32(0)


process.demo = cms.EDAnalyzer('DemoAnalyzer')



# ***********************************************************
# output file name                                          *
# default is DoubleMuParked2012C_10000_Higgs.root           *
# ***********************************************************
#process.TFileService = cms.Service("TFileService",
#       fileName = cms.string('Practice.root')
#                                   )

process.p = cms.Path(process.demo)
