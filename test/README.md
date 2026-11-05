Test
====

    cmsRun runDump.py  inputFiles=file:/tmp/TEST_UncalibRecHitsOut.root   outputFile=test.root   maxEvents=20

    cmsRun runDump.py  inputFiles=file:/tmp/TEST_UncalibRecHitsOut.root   outputFile=test.root



Phase-2
====


    dasgoclient --query "file dataset=/RelValSingleEFlatPt2To100/treis-RelValSingleEFlatPt2To100_nopu_ecal_612wf_step2_2025ps-afa0981336aeeaf41d23e9ebb0443183/USER instance=prod/phys03"

    /store/user/treis/RelValSingleEFlatPt2To100/RelValSingleEFlatPt2To100_nopu_ecal_612wf_step2_2025ps/261009_085711/0000/step2_3.root
    /store/user/treis/RelValSingleEFlatPt2To100/RelValSingleEFlatPt2To100_nopu_ecal_612wf_step2_2025ps/261009_085711/0000/step2_1.root
    /store/user/treis/RelValSingleEFlatPt2To100/RelValSingleEFlatPt2To100_nopu_ecal_612wf_step2_2025ps/261009_085711/0000/step2_2.root

    edmDumpEventContent /store/user/treis/RelValSingleEFlatPt2To100/RelValSingleEFlatPt2To100_nopu_ecal_612wf_step2_2025ps/261009_085711/0000/step2_2.root > dump.txt

    edmDumpEventContent /store/user/treis/RelValSingleEFlatPt2To100/RelValSingleEFlatPt2To100_nopu_ecal_612wf_step2_2025ps/261009_085711/0000/step2_3.root > dump_3.txt



    cmsRun runSamplesDump.py inputFiles=file:/store/user/treis/RelValSingleEFlatPt2To100/RelValSingleEFlatPt2To100_nopu_ecal_612wf_step2_2025ps/261009_085711/0000/step2_3.root outputFile=/tmp/amassiro/samples.root

    cmsRun runSamplesDump.py inputFiles=file:/eos/cms/store/user/treis/RelValSingleEFlatPt2To100/RelValSingleEFlatPt2To100_nopu_ecal_612wf_step2_2025ps/261009_085711/0000/step2_3.root outputFile=/tmp/amassiro/samples.root


Plot the pulses
====

Pedestal = (ADC[0] + ADC[1]) / 2, subtracted to all the samples. Pulses normalized to ADC[5] - pedestal, cut ADC[5] - pedestal > threshold.

One crystal (ieta, iphi), one event or mean over all events (iEvent = -1):

    root -l 'drawPulse.cxx("samples.root", -54, 119, 3)'
    root -l 'drawPulse.cxx("samples.root", -54, 119, -1, 50)'

One ieta, all iphi:

    root -l 'drawPulseEta.cxx("samples.root", -54, 50)'

All ieta (-85 ... 85), plots in pulses_output/ (.png and .root), summary plot and mean pulses in pulses_output/pulses.txt (ieta, ADC[0], ... , ADC[15]):

    root -l -b -q 'drawPulseAllEta.cxx+("samples.root", 50)'

Weights from the mean pulses (ieta = 0 means all ieta):

    root -l -q 'LocalCreateWeightsFromPulses.cxx("pulses_output/pulses.txt", 0, "weights_pulses_all.txt")'

Apply the weights:

    root -l 'ApplyWeights.cxx("samples.root", "weights_pulses_all.txt", 50, -54, 119)'

