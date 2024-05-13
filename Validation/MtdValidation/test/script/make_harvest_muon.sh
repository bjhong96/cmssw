#!/bin/sh

dirname=$1

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

linenum_tot=$(awk "END { print NR }" ${dirname}/flist.txt)
linenum=0


cat << 'EOF_1' > ${dirname}/mtdHarvesting_muon_cfg_${linenum_tot}.py
import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils


from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9
process = cms.Process('mtdHarvesting',Phase2C17I13M9)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.EDMtoMEAtRunEnd_cff')       # This line is important
process.load('SimGeneral.MixingModule.mixNoPU_cfi')

process.load("Configuration.Geometry.GeometryExtended2026D98Reco_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#mylist = FileUtils.loadListFromFile('filenames_update.txt') # input file with file names from grid grid control output

process.MessageLogger.cerr.FwkReport  = cms.untracked.PSet(
    reportEvery = cms.untracked.int32(-1),
)

# Input source
process.source = cms.Source("DQMRootSource",
    fileNames = cms.untracked.vstring(
EOF_1
while read line; do
  ((linenum+=1))
  if [[ "$linenum" == "$linenum_tot" ]]; then
    echo "\"file:output/step3_inDQM_${p_or_np}_${linenum_tot}_${linenum}.root\"" >> ${dirname}/mtdHarvesting_muon_cfg_${linenum_tot}.py
  else
    echo "\"file:output/step3_inDQM_${p_or_np}_${linenum_tot}_${linenum}.root\"," >> ${dirname}/mtdHarvesting_muon_cfg_${linenum_tot}.py
  fi
done < ${dirname}/flist.txt
cat << "EOF_2" >> ${dirname}/mtdHarvesting_muon_cfg_${linenum_tot}.py
    )
)

# input source, when using grid control - from file with filename list
# process.source = cms.Source("DQMRootSource",
#     fileNames = cms.untracked.vstring(*mylist)
# )

# Path and EndPath definitions

process.edmtome_step = cms.Path(process.EDMtoME)
process.dqmsave_step = cms.Path(process.DQMSaver)

# --- PostProcessing

process.load("Validation.MtdValidation.btlSimHitsPostProcessor_cfi")
process.load("Validation.MtdValidation.btlLocalRecoPostProcessor_cfi")
process.load("Validation.MtdValidation.MtdTracksPostProcessor_cfi")
#process.load("Validation.MtdValidation.MtdEleIsoPostProcessor_cfi")   # Normunds
process.load("Validation.MtdValidation.MtdMuonIsoPostProcessor_cfi")
process.load("Validation.MtdValidation.Primary4DVertexPostProcessor_cfi")

#process.harvesting = cms.Sequence(process.btlSimHitsPostProcessor + process.btlLocalRecoPostProcessor + process.MtdTracksPostProcessor + process.Primary4DVertexPostProcessor)
#process.harvesting = cms.Sequence(process.btlSimHitsPostProcessor + process.btlLocalRecoPostProcessor + process.MtdTracksPostProcessor + process.MtdEleIsoPostProcessor + process.Primary4DVertexPostProcessor)   # Normunds
#process.harvesting = cms.Sequence(process.btlSimHitsPostProcessor + process.btlLocalRecoPostProcessor + process.MtdTracksPostProcessor + process.MtdEleIsoPostProcessor + process.MtdMuonIsoPostProcessor + process.Primary4DVertexPostProcessor)
process.harvesting = cms.Sequence(process.btlSimHitsPostProcessor + process.btlLocalRecoPostProcessor + process.MtdTracksPostProcessor + process.MtdMuonIsoPostProcessor + process.Primary4DVertexPostProcessor)

process.p = cms.Path( process.harvesting )

process.schedule = cms.Schedule( process.edmtome_step , process.p , process.dqmsave_step )
EOF_2
