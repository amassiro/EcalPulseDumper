# EcalPulseDumper

Where:

    /afs/cern.ch/user/a/amassiro/work/ECAL/GPU/onGPU/24Jun2021/CMSSW_11_3_0_pre5/src/ECALValidation

To install:

    mkdir ECALValidation
    cd ECALValidation/
    git clone git@github.com:amassiro/EcalPulseDumper.git

    cmsenv
    scramv1 b -j 20


    
    

Phase-2
====

SamplesTreeProducer: dump only the ECAL phase-2 barrel digi samples (EBDigiCollectionPh2):
16 samples (6.25 ns), raw ADC (12 bits) and gain id (0 = x10, 1 = x1), ieta, iphi.
Optionally amplitude and time of the CMSSW uncalibrated rechits (ecalUncalibRecHitPhase2, weights).

To install: as above, in a recent release with phase-2 ECAL digis (e.g. CMSSW_15_1_X, used for the TB2025 pulse shape tag)

Run:

    cd test
    cmsRun runSamplesDump.py inputFiles=file:sample.root outputFile=samples.root maxEvents=100
    cmsRun runSamplesDump.py inputFiles=file:sample.root outputFile=samples.root digiTag=simEcalUnsuppressedDigis minAmplitude=10
    cmsRun runSamplesDump.py inputFiles=file:sample.root outputFile=samples.root rechitTag=ecalUncalibRecHitPhase2:EcalUncalibRecHitsEB

Where:

    digiTag      = simEcalUnsuppressedDigis (mixData with premixing), check with edmDumpEventContent
    rechitTag    = only if the sample has phase-2 uncalibrated rechits (era with phase2_ecal_devel)
    minAmplitude = keep digis with max(ADC) - ADC[0] > minAmplitude, or with gain x1 (0 = keep all)

NB: to simulate with the TB2025 pulse shape add to cmsDriver:

    --custom_conditions EcalSimPulseShapes_TB2025_v1,EcalSimPulseShapeRcd,frontier://FrontierProd/CMS_CONDITIONS

Apply the weights (from EcalWeights/LocalCreateWeightsPhase2.cxx):

    root -l 'ApplyWeights.cxx("samples.root", "weights_TB2025_v1.txt")'

