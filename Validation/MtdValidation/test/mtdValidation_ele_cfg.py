
# imports our CMS-specific Python classes and functions
import FWCore.ParameterSet.Config as cms

# Eras control the detector scenario
from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9
# the Process
process = cms.Process('mtdValidation',Phase2C17I13M9)

process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')

# GEANT simulation needs geometry and magnetic field
process.load("Configuration.Geometry.GeometryExtended2026D95Reco_cff")
process.load('Configuration.StandardSequences.MagneticField_cff')
# geometry and magnetic are stored in database
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
# global tag needs to match the Era
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T21', '')

process.load('RecoLocalFastTime.FTLClusterizer.MTDCPEESProducer_cfi')
process.load("Configuration.StandardSequences.Reconstruction_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#Setup FWK for multithreaded    FWK: Framework
process.options.numberOfThreads = 4
process.options.numberOfStreams = 0
process.options.numberOfConcurrentLuminosityBlocks = 0
process.options.eventSetup.numberOfConcurrentIOVs = 1

process.MessageLogger.cerr.FwkReport  = cms.untracked.PSet(
    reportEvery = cms.untracked.int32(100),
)

# input files from correct CMSSW version!
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#        'file:step3.root'
#	'/store/relval/CMSSW_13_1_0_pre3/RelValZEE_14/GEN-SIM-RECO/PU_131X_mcRun4_realistic_v2_2026D95PU200-v1/00000/00ac807f-b506-4e65-8f14-49f9f9e9dbf5.root'
#	'file:/eos/cms/store/relval/CMSSW_13_2_0_pre1/RelValZMM_14/GEN-SIM-RECO/PU_131X_mcRun4_realistic_v5_2026D95PU200-v1/2590000/0d9804a2-bdf6-45a6-8a97-8440b3dc88fa.root'
#	'file:/eos/cms/store/relval/CMSSW_13_2_0_pre1/RelValZEE_14/GEN-SIM-RECO/PU_131X_mcRun4_realistic_v5_2026D95PU200-v1/2590000/0ed9d685-abfc-4b43-8b89-26099a64f09a.root'
	'file:/eos/cms/store/relval/CMSSW_13_2_0_pre1/RelValTTbar_14TeV/GEN-SIM-RECO/PU_131X_mcRun4_realistic_v5_2026D95PU200-v1/2590000/237052f0-eff3-431f-81a2-9e34823953d7.root'
    )
)

process.mix.digitizers = cms.PSet()
for a in process.aliases: delattr(process, a)

# --- BTL Validation
process.load("Validation.MtdValidation.btlSimHitsValid_cfi")
process.load("Validation.MtdValidation.btlDigiHitsValid_cfi")
process.load("Validation.MtdValidation.btlLocalRecoValid_cfi")
btlValidation = cms.Sequence(process.btlSimHitsValid + process.btlDigiHitsValid + process.btlLocalRecoValid)

# --- ETL Validation
process.load("Validation.MtdValidation.etlSimHitsValid_cfi")
process.load("Validation.MtdValidation.etlDigiHitsValid_cfi")
process.load("Validation.MtdValidation.etlLocalRecoValid_cfi")
etlValidation = cms.Sequence(process.etlSimHitsValid + process.etlDigiHitsValid + process.etlLocalRecoValid)

# --- Global Validation
process.load("Validation.MtdValidation.mtdTracksValid_cfi")
process.load("Validation.MtdValidation.mtdEleIsoValid_cfi")   #Normunds
process.load("Validation.MtdValidation.vertices4DValid_cfi")

# process.btlDigiHitsValid.optionalPlots = True
# process.etlDigiHitsValid.optionalPlots = True
# process.btlLocalRecoValid.optionalPlots = True
# process.etlLocalRecoValid.optionalPlots = True
# process.mtdTracksValid.optionalPlots = True
# process.vertices4DValid.optionalPlots = True

#process.validation = cms.Sequence(btlValidation + etlValidation + process.mtdTracksValid + process.vertices4DValid)
process.validation = cms.Sequence(btlValidation + etlValidation + process.mtdTracksValid + process.mtdEleIsoValid + process.vertices4DValid)   #Normunds

process.DQMoutput = cms.OutputModule("DQMRootOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('DQMIO'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:step3_inDQM.root'),
    outputCommands = process.DQMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

process.p = cms.Path( process.mix + process.mtdTrackingRecHits + process.validation )
process.endjob_step = cms.EndPath(process.endOfProcess)
process.DQMoutput_step = cms.EndPath( process.DQMoutput )

process.schedule = cms.Schedule( process.p , process.endjob_step , process.DQMoutput_step )
