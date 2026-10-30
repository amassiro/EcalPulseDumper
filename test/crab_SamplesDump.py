# CRAB config to run SamplesTreeProducer on the phase-2 RelVal (step2 = DIGI, private production in phys03)
#
#   crab submit -c crab_SamplesDump.py
#   crab status -d crab_projects/crab_SingleEFlatPt2To100_2025ps

from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.requestName  = 'SingleEFlatPt2To100_2025ps'
config.General.workArea     = 'crab_projects'
config.General.transferLogs = True

config.section_('JobType')
config.JobType.pluginName   = 'Analysis'
config.JobType.psetName     = 'runSamplesDump.py'
config.JobType.pyCfgParams  = ['outputFile=samples.root', 'runWeights=True']
config.JobType.outputFiles  = ['samples.root']
config.JobType.maxMemoryMB  = 2500
config.JobType.allowUndistributedCMSSW = True

config.section_('Data')
config.Data.inputDataset  = '/RelValSingleEFlatPt2To100/treis-RelValSingleEFlatPt2To100_nopu_ecal_612wf_step2_2025ps-afa0981336aeeaf41d23e9ebb0443183/USER'
config.Data.inputDBS      = 'phys03'
config.Data.splitting     = 'FileBased'
config.Data.unitsPerJob   = 5
config.Data.outLFNDirBase = '/store/user/amassiro/ECAL/Phase2'
config.Data.publication   = False

config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'
