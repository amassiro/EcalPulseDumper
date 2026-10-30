// -*- C++ -*-
//
// Package:    EcalPulseDumper
// Class:      SamplesTreeProducer
//
/**\class SamplesTreeProducer SamplesTreeProducer.cc EcalPulseDumper/plugins/SamplesTreeProducer.cc

 Description: dump only the digi samples of ECAL phase-2 barrel (EBDigiCollectionPh2)

 Implementation:
   - 16 samples per digi, 160 MHz (6.25 ns), see DataFormats/EcalDigi/interface/EcalConstants.h
   - for each sample: raw ADC (12 bits) and gain id (0 = x10, 1 = x1)
   - only digis present in the collection are stored (variable size arrays),
     optionally only if max(ADC) - ADC[0] > minAmplitude
   - optionally (EcalUncalibRecHitsEBCollection not empty) amplitude and time
     of the CMSSW uncalibrated rechit of the same crystal, e.g. ecalUncalibRecHitPhase2 (weights)
*/
//
// Original Author:  Andrea Massironi
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"


// ECAL specific

#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalDigi/interface/EcalDataFrame_Ph2.h"
#include "DataFormats/EcalDigi/interface/EcalLiteDTUSample.h"
#include "DataFormats/EcalDigi/interface/EcalConstants.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalUncalibratedRecHit.h"


#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "TTree.h"

//
// class declaration
//

class SamplesTreeProducer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit SamplesTreeProducer(const edm::ParameterSet& );
  ~SamplesTreeProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void analyze(const edm::Event&, const edm::EventSetup&) override;

  static constexpr int kNSamples  = ecalPh2::sampleSize;     // 16
  static constexpr int kNChannels = ecalPh2::kEBChannels;    // 61200

  // ----------member data ---------------------------

  const edm::EDGetTokenT<EBDigiCollectionPh2> _token_ebdigi;
  const bool _doRecHits;
  edm::EDGetTokenT<EcalUncalibratedRecHitCollection> _token_ebrechits;
  const int _minAmplitude;

  TTree *_outTree;

  UInt_t    _run;
  UShort_t  _lumi;
  UShort_t  _bx;
  ULong64_t _event;

  int _nDigis;       // number of stored digis
  int _nSamplesTot;  // _nDigis * kNSamples

  int      _ieta[kNChannels];
  int      _iphi[kNChannels];
  UShort_t _adc[kNChannels * kNSamples];
  UChar_t  _gain[kNChannels * kNSamples];

  float    _amplitude_cmssw[kNChannels];   // -999 if no rechit
  float    _time_cmssw[kNChannels];

};


//
// constructors and destructor
//
SamplesTreeProducer::SamplesTreeProducer(const edm::ParameterSet& iConfig) :
  _token_ebdigi(consumes<EBDigiCollectionPh2>(iConfig.getParameter<edm::InputTag>("EBDigiCollection"))),
  _doRecHits(!iConfig.getParameter<edm::InputTag>("EcalUncalibRecHitsEBCollection").label().empty()),
  _minAmplitude(iConfig.getParameter<int>("minAmplitude"))
{
  if (_doRecHits) _token_ebrechits = consumes<EcalUncalibratedRecHitCollection>(iConfig.getParameter<edm::InputTag>("EcalUncalibRecHitsEBCollection"));

  usesResource("TFileService");
  edm::Service<TFileService> fs;

  _outTree = fs->make<TTree>("tree","tree");

  _outTree->Branch("run",           &_run,           "run/i");
  _outTree->Branch("lumi",          &_lumi,          "lumi/s");
  _outTree->Branch("bx",            &_bx,            "bx/s");
  _outTree->Branch("event",         &_event,         "event/l");

  _outTree->Branch("nDigis",        &_nDigis,        "nDigis/I");
  _outTree->Branch("nSamplesTot",   &_nSamplesTot,   "nSamplesTot/I");

  _outTree->Branch("ieta",          _ieta,           "ieta[nDigis]/I");
  _outTree->Branch("iphi",          _iphi,           "iphi[nDigis]/I");
  _outTree->Branch("adc",           _adc,            "adc[nSamplesTot]/s");   // adc [iDigi*16 + iSample]
  _outTree->Branch("gain",          _gain,           "gain[nSamplesTot]/b");  // 0 = x10, 1 = x1

  if (_doRecHits) {
    _outTree->Branch("amplitude_cmssw", _amplitude_cmssw, "amplitude_cmssw[nDigis]/F");
    _outTree->Branch("time_cmssw",      _time_cmssw,      "time_cmssw[nDigis]/F");
  }

}


//
// member functions
//

// ------------ method called for each event  ------------
void
SamplesTreeProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  _run   = iEvent.eventAuxiliary().run();
  _lumi  = iEvent.eventAuxiliary().luminosityBlock();
  _bx    = iEvent.eventAuxiliary().bunchCrossing();
  _event = iEvent.eventAuxiliary().event();

  const EBDigiCollectionPh2& ebdigis = iEvent.get(_token_ebdigi);

  const EcalUncalibratedRecHitCollection* ebrechits = nullptr;
  if (_doRecHits) ebrechits = &iEvent.get(_token_ebrechits);

  _nDigis = 0;

  //---- dump digis
  for (auto itdigi = ebdigis.begin(); itdigi != ebdigis.end(); itdigi++) {

    EcalDataFrame_Ph2 dataFrame(*itdigi);
    if (dataFrame.size() != kNSamples) {
      edm::LogWarning("SamplesTreeProducer") << " digi with " << dataFrame.size() << " samples, expected " << kNSamples << ": skipped";
      continue;
    }

    //---- optional selection: max(ADC) - ADC[0] > minAmplitude, only for gain x10 samples
    if (_minAmplitude > 0) {
      int first_adc = dataFrame[0].adc();
      int max_adc = first_adc;
      bool has_gain1 = false;
      for (int iSample = 0; iSample < kNSamples; iSample++) {
        if (dataFrame[iSample].gainId() == int(ecalPh2::gainId1)) has_gain1 = true;
        if (dataFrame[iSample].adc() > max_adc) max_adc = dataFrame[iSample].adc();
      }
      if (!has_gain1 && (max_adc - first_adc) <= _minAmplitude) continue;
    }

    EBDetId id(dataFrame.id());
    _ieta[_nDigis] = id.ieta();
    _iphi[_nDigis] = id.iphi();

    if (_doRecHits) {
      _amplitude_cmssw[_nDigis] = -999;
      _time_cmssw[_nDigis] = -999;
      auto itrechit = ebrechits->find(id);
      if (itrechit != ebrechits->end()) {
        _amplitude_cmssw[_nDigis] = itrechit->amplitude();
        _time_cmssw[_nDigis] = itrechit->jitter();
      }
    }

    for (int iSample = 0; iSample < kNSamples; iSample++) {
      EcalLiteDTUSample sample = dataFrame[iSample];
      _adc [_nDigis * kNSamples + iSample] = sample.adc();
      _gain[_nDigis * kNSamples + iSample] = sample.gainId();
    }

    _nDigis++;
  }

  _nSamplesTot = _nDigis * kNSamples;

  _outTree->Fill();

}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SamplesTreeProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("EBDigiCollection", edm::InputTag("simEcalUnsuppressedDigis", ""));
  desc.add<edm::InputTag>("EcalUncalibRecHitsEBCollection", edm::InputTag(""));  // e.g. ecalUncalibRecHitPhase2:EcalUncalibRecHitsEB, empty = not dumped
  desc.add<int>("minAmplitude", 0);  // in ADC counts, 0 = keep all digis
  descriptions.add("SamplesTreeProducer", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SamplesTreeProducer);

