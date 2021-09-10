# test reco and dump into a tree

import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('analysis')
options.parseArguments()

process = cms.Process('ECALNoise')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(10000)

process.TFileService = cms.Service("TFileService",
     fileName = cms.string(options.outputFile)
)

process.options = cms.untracked.PSet(
#    SkipEvent = cms.untracked.vstring('ProductNotFound'),
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('reco nevts:1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag


# 2021 ...
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2021_realistic', '')

process.PulseTreeProducer = cms.EDAnalyzer('PulseTreeProducer',
                            EcalUncalibRecHitsEBCollection = cms.InputTag("ecalMultiFitUncalibRecHit",  "EcalUncalibRecHitsEB"),
                            EcalUncalibRecHitsEECollection = cms.InputTag("ecalMultiFitUncalibRecHit",  "EcalUncalibRecHitsEE"),

                            EcalUncalibRecHitsEBCollectionSecond = cms.InputTag("ecalUncalibRecHitConvertGPU2CPUFormat", "EcalUncalibRecHitsEBfromGPU"),
                            EcalUncalibRecHitsEECollectionSecond = cms.InputTag("ecalUncalibRecHitConvertGPU2CPUFormat", "EcalUncalibRecHitsEEfromGPU"),

                            #EcalRecHitsEBCollection = cms.InputTag("ecalRecHit",  "EcalRecHitsEB"),
                            #EcalRecHitsEECollection = cms.InputTag("ecalRecHit",  "EcalRecHitsEE"),

                            EBDigiCollection = cms.InputTag("ecalDigis",  "ebDigis"),
                            EEDigiCollection = cms.InputTag("ecalDigis",  "eeDigis"),

                           )

#edm::SortedCollection<EcalUncalibratedRecHit,edm::StrictWeakOrdering<EcalUncalibratedRecHit> >    "ecalMultiFitUncalibRecHit"   "EcalUncalibRecHitsEB"   "RECO"    
#EBDigiCollection                      "ecalDigis"                 "ebDigis"         "RECO"    


process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(50)

process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(options.inputFiles),
                                secondaryFileNames = cms.untracked.vstring()
                                )


process.PulseTreeProducer_step = cms.Path(process.PulseTreeProducer)
process.endjob_step = cms.EndPath(process.endOfProcess)


process.schedule = cms.Schedule(
                                process.PulseTreeProducer_step, 
                                process.endjob_step
                                )


