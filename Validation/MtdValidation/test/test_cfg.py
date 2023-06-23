import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9
process = cms.Process('mtdHarvesting',Phase2C17I13M9)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.EDMtoMEAtRunEnd_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load("Configuration.Geometry.GeometryExtended2026D95Reco_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.MessageLogger.cerr.FwkReport = cms.untracked.PSet(
	reportEvery = cms.untracked.int32(-1),
)

process.source = cms.Source("DQMRootSource",
			fileNames = cms.untracked.vstring("file:step3_inDQM.root")
			)

process.edmtome_step = cms.Path(process.EDMtoME)
process.dqmsave_step = cms.Path(process.DQMSaver)  # This line makes an output file

process.load("Validation.MtdValidation.btlSimHitsPostProcessor_cfi")

process.harvesting = cms.Sequence(process.btlSimHitsPostProcessor)

process.p = cms.Path( process.harvesting )    # This line cout "hi"
