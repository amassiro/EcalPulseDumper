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


Test on all samples
====

    cmsRun runSamplesDump.py inputFiles=file:/eos/cms/store/user/treis/RelValSingleEFlatPt2To100/RelValSingleEFlatPt2To100_nopu_ecal_612wf_step2_2025ps/261009_085711/0000/step2_1.root outputFile=/tmp/amassiro/samples_1.root
    cmsRun runSamplesDump.py inputFiles=file:/eos/cms/store/user/treis/RelValSingleEFlatPt2To100/RelValSingleEFlatPt2To100_nopu_ecal_612wf_step2_2025ps/261009_085711/0000/step2_2.root outputFile=/tmp/amassiro/samples_2.root
    cmsRun runSamplesDump.py inputFiles=file:/eos/cms/store/user/treis/RelValSingleEFlatPt2To100/RelValSingleEFlatPt2To100_nopu_ecal_612wf_step2_2025ps/261009_085711/0000/step2_3.root outputFile=/tmp/amassiro/samples_3.root

Merge:

    hadd samples_123.root /tmp/amassiro/samples_1.root /tmp/amassiro/samples_2.root /tmp/amassiro/samples_3.root


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

All crystals and all events together, plots in pulses_output/pulse_average and pulse_2d_average (.png and .root):

    root -l 'drawPulseAverage.cxx+("samples.root", 50)'

Weights from this mean pulse (pulses_output/pulses_average.txt), amplitude + pedestal only, no time, no noise correlation,
i.e. w_i = (f_i - sum_j f_j / N) / (sum_j f_j^2 - (sum_j f_j)^2 / N):

    root -l -q 'LocalCreateWeightsFromPulses.cxx("pulses_output/pulses_average.txt", 0, "weights_pulses_average_simple.txt", true, false, false)'


Phase-2 weights for CMSSW
====

Amplitude and time weights in RecoLocalCalo/EcalRecProducers/plugins/EcalUncalibRecHitPhase2WeightsProducer.cc
(branch amassiro-phase2-weights of github.com/amassiro/cmssw).

Where:

    /home/amassiro/Cern/Claude/ECAL/data/samples.root = hadd of samples_1.root, samples_2.root, samples_3.root (9000 events)

Mean pulse over all crystals and all events, ADC[5] - pedestal > 50 ADC (23845 pulses):

    root -l -b -q 'drawPulseAverage.cxx+("../../data/samples.root", 50, "pulses_output_data")'

Weights: least squares fit of amplitude A, A * dt and pedestal P (doFitBaseline = true, doFitTime = true, useNoiseCorrelation = false),
W = (C^T C)^-1 C^T, C columns = mean pulse f, derivative df/dt (cubic spline), 1 (CMS NOTE 2006/037, eq. A.5 - A.6).
The ieta of the mean pulse is renamed from 0 to 99, because ieta = 0 means "all ieta" and does not print the time weights:

    sed 's/^0 /99 /' pulses_output_data/pulses_average.txt > pulses_output_data/pulses_average_99.txt
    root -l -q 'LocalCreateWeightsFromPulses.cxx("pulses_output_data/pulses_average_99.txt", 99, "weights_pulses_average_data.txt", true, true, false)'

ampWeights and timeWeights are printed on screen as cms.vdouble(...), to be copied in the producer.
Amplitude: sum_i w_i S_i = A. Time: sum_i t_i S_i = A * dt [ADC * ns], dt > 0 means earlier pulse.

NB: the values in the branch amassiro-phase2-weights are from samples.root = step2_3 only (4542 pulses):

    root -l -b -q 'drawPulseAverage.cxx+("samples.root", 50)'
    sed 's/^0 /99 /' pulses_output/pulses_average.txt > pulses_output/pulses_average_99.txt
    root -l -q 'LocalCreateWeightsFromPulses.cxx("pulses_output/pulses_average_99.txt", 99, "weights_pulses_average_time.txt", true, true, false)'

Weights from the mean pulses (ieta = 0 means all ieta):

    root -l -q 'LocalCreateWeightsFromPulses.cxx("pulses_output/pulses.txt", 0, "weights_pulses_all.txt")'

Weights per crystal, from samples.root (threshold on ADC[5] - pedestal, only crystals with at least 3 pulses, their number is printed), output "ieta iphi w[0] ... w[15]":

    root -l -q 'LocalCreateWeightsFromPulsesAll.cxx("samples.root", 50, "weights_pulses_crystals.txt")'

Average of the weights over the ieta rings, in average_weights_pulses_all.txt, and distributions of w[0] ... w[15] in pulses_output/weights_distributions.png (.root):

    root -l 'averageWeights.cxx("weights_pulses_all.txt")'

Same, with the weights per crystal, in average_weights_pulses_crystals.txt and pulses_output/weights_distributions_crystals.png (.root):

    root -l 'averageWeightsCrystals.cxx("weights_pulses_crystals.txt")'

Apply the weights:

    root -l 'ApplyWeights.cxx("samples.root", "weights_pulses_all.txt", 50, -54, 119)'

With the unique weights for all crystals (file with one line only, e.g. from the mean pulse of all crystals and all events):

    root -l 'ApplyWeights.cxx("../../data/samples.root", "weights_pulses_average_data.txt", 50, -54, 119)'

