# dump only the ECAL phase-2 barrel digi samples into a tree

import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('analysis')
options.register ('digiTag',
                  'simEcalUnsuppressedDigis',
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  "EBDigiCollectionPh2 input tag")
options.register ('rechitTag',
                  '',
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  "EcalUncalibRecHitCollection EB input tag, e.g. ecalUncalibRecHitPhase2:EcalUncalibRecHitsEB (empty = not dumped)")
options.register ('runWeights',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.bool,
                  "rerun ecalUncalibRecHitPhase2 (CMSSW weights) on the digis and dump its amplitude and time")
options.register ('minAmplitude',
                  0,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.int,
                  "keep digis with max(ADC) - ADC[0] > minAmplitude (0 = keep all)")
options.parseArguments()

process = cms.Process('ECALSamples')

process.load('FWCore.MessageService.MessageLogger_cfi')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(50)

process.TFileService = cms.Service("TFileService",
     fileName = cms.string(options.outputFile)
)

process.SamplesTreeProducer = cms.EDAnalyzer('SamplesTreeProducer',
                            EBDigiCollection = cms.InputTag(options.digiTag),
                            EcalUncalibRecHitsEBCollection = cms.InputTag(options.rechitTag),
                            minAmplitude = cms.int32(options.minAmplitude),
                           )

process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(options.inputFiles),
                                secondaryFileNames = cms.untracked.vstring()
                                )

process.SamplesTreeProducer_step = cms.Path(process.SamplesTreeProducer)

#
# rerun the CMSSW phase-2 weights on the digis (no EventSetup needed, weights hardcoded in the producer)
#
if options.runWeights:
  process.load('RecoLocalCalo.EcalRecProducers.ecalUncalibRecHitPhase2_cfi')
  process.ecalUncalibRecHitPhase2.BarrelDigis = cms.InputTag(options.digiTag)
  process.SamplesTreeProducer.EcalUncalibRecHitsEBCollection = cms.InputTag('ecalUncalibRecHitPhase2', 'EcalUncalibRecHitsEB')
  process.SamplesTreeProducer_step = cms.Path(process.ecalUncalibRecHitPhase2 * process.SamplesTreeProducer)

process.schedule = cms.Schedule(
                                process.SamplesTreeProducer_step,
                                )
