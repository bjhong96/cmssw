#!/bin/sh

path=$1
dirname=$2
#dirname="240227_PU200_prompt_muon_test"
mkdir -p ${dirname}
mkdir -p ${dirname}/output
mkdir -p ${dirname}/log

if [[ "$dirname" == *"_prompt"* ]]; then
  echo "samples for prompt muon"
  p_or_np="prompt"
elif [[ "$dirname" == *"_nonprompt"* ]]; then
  echo "samples for nonprompt muon"
  p_or_np="nonprompt"
else
  echo "dirname should include \"prompt\" or \"nonprompt\""
  break
fi

if [[ "$dirname" == *"_PU200"* ]]; then
  echo "samples for PU200 scenario"
  pu_or_nopu="PU200"
elif [[ "$dirname" == *"_noPU"* ]]; then
  echo "samples for noPU scenario"
  pu_or_nopu="noPU"
else
  echo "dirname should include \"PU200\" or \"noPU\""
  break
fi

find ${path} -type f > ${dirname}/flist.txt
sleep 5s
sed -i 's/^/"file:/' ${dirname}/flist.txt
sed -i 's/root/root"/' ${dirname}/flist.txt

linenum_tot=$(awk "END { print NR }" ${dirname}/flist.txt)
linenum=0

while read line; do
  ((linenum+=1))
  echo $linenum

  cat << 'EOF_1' > ${dirname}/mtdValidation_muon_cfg_${p_or_np}_${linenum_tot}_${linenum}.py
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
process.load("Configuration.Geometry.GeometryExtended2026D98Reco_cff")
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
#process.options.numberOfThreads = 4
#process.options.numberOfStreams = 0
process.options.numberOfConcurrentLuminosityBlocks = 0
process.options.eventSetup.numberOfConcurrentIOVs = 1

process.MessageLogger.cerr.FwkReport  = cms.untracked.PSet(
    reportEvery = cms.untracked.int32(100),
)

# input files from correct CMSSW version!
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
EOF_1
  echo $line >> ${dirname}/mtdValidation_muon_cfg_${p_or_np}_${linenum_tot}_${linenum}.py
  cat << "EOF_2" >> ${dirname}/mtdValidation_muon_cfg_${p_or_np}_${linenum_tot}_${linenum}.py
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
#process.load("Validation.MtdValidation.mtdEleIsoValid_cfi)    # Electron
process.load("Validation.MtdValidation.mtdMuonIsoValid_cfi")
process.load("Validation.MtdValidation.vertices4DValid_cfi")

# process.btlDigiHitsValid.optionalPlots = True
# process.etlDigiHitsValid.optionalPlots = True
# process.btlLocalRecoValid.optionalPlots = True
# process.etlLocalRecoValid.optionalPlots = True
# process.mtdTracksValid.optionalPlots = True
# process.vertices4DValid.optionalPlots = True

#process.validation = cms.Sequence(btlValidation + etlValidation + process.mtdTracksValid + process.vertices4DValid)
#process.validation = cms.Sequence(btlValidation + etlValidation + process.mtdTracksValid + process.mtdEleIsoValid + process.mtdMuonIsoValid + process.vertices4DValid)
process.validation = cms.Sequence(btlValidation + etlValidation + process.mtdTracksValid + process.mtdMuonIsoValid + process.vertices4DValid)

process.DQMoutput = cms.OutputModule("DQMRootOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('DQMIO'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string(
EOF_2
  echo "\"file:${PWD}/${dirname}/output/step3_inDQM_${p_or_np}_${linenum_tot}_${linenum_tot}_${linenum}.root\"" >> ${dirname}/mtdValidation_muon_cfg_${p_or_np}_${linenum_tot}_${linenum}.py
  cat << "EOF_3" >> ${dirname}/mtdValidation_muon_cfg_${p_or_np}_${linenum_tot}_${linenum}.py
    ),
    #fileName = cms.untracked.string('file:step3_inDQM_sig_48.root'),
    outputCommands = process.DQMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)
# Get Ntuples
#process.TFileService = cms.Service("TFileService", fileName=cms.string('/afs/cern.ch/user/b/byhong/CMSSW_13_3_0_pre1/src/Validation/MtdValidation/test/231123_PU200_prompt/output/ntuple_prompt_180.root'))

process.p = cms.Path( process.mix + process.mtdTrackingRecHits + process.validation )
process.endjob_step = cms.EndPath(process.endOfProcess)
process.DQMoutput_step = cms.EndPath( process.DQMoutput )

process.schedule = cms.Schedule( process.p , process.endjob_step , process.DQMoutput_step )

EOF_3

done < ${dirname}/flist.txt
