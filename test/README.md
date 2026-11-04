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
