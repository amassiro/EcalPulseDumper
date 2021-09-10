// -*- C++ -*-
//
// Package:    EcalZee
// Class:      PulseTreeProducer
// 
/**\class PulseTreeProducer PulseTreeProducer.cc EcalZee/plugins/PulseTreeProducer.cc
 * 
 Description: [one line class summary]
 
 Implementation:
 [Notes on implementation]
 */
//
// Original Author:  Andrea Massironi
//         Created:  Thu, 1 July 2021 10:09:05 GMT
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


#include "FWCore/Framework/interface/EventSetup.h"




// ECAL specific

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalUncalibratedRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"

#include "CondFormats/EcalObjects/interface/EcalPedestals.h"
#include "CondFormats/DataRecord/interface/EcalPedestalsRcd.h"




#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"


#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "TTree.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class PulseTreeProducer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit PulseTreeProducer(const edm::ParameterSet& ); // z, edm::ConsumesCollector& );
  ~PulseTreeProducer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  // ----------member data ---------------------------
  
  edm::EDGetTokenT<EcalUncalibratedRecHitCollection> _token_ebrechits;
  edm::EDGetTokenT<EcalUncalibratedRecHitCollection> _token_eerechits;
  
  edm::EDGetTokenT<EcalUncalibratedRecHitCollection> _token_second_ebrechits;
  edm::EDGetTokenT<EcalUncalibratedRecHitCollection> _token_second_eerechits;
  
  edm::EDGetTokenT<EBDigiCollection> _token_ebdigi;
  edm::EDGetTokenT<EEDigiCollection> _token_eedigi;
  
  
//   edm::ESHandle<EcalPedestals> _peds;
//   edm::ESGetToken<EcalPedestals, EcalPedestalsRcd> _pedsToken;
  
  const EcalPedestals* _peds;
  
  
  TTree *_outTree;
  
  UInt_t _run;
  UShort_t _lumi;
  UShort_t _bx;
  UShort_t _event;      
  
  
  float _amplitude_EB[61200];
  float _amplitude_second_EB[61200];
  float _chi2_EB[61200];
  float _chi2_second_EB[61200];
  float _amplitudeError_EB[61200];
  float _amplitudeError_second_EB[61200];
  int   _ieta[61200];
  int   _iphi[61200];
  float _digi_ped_subtracted_EB[61200*10];
  
  
  float _amplitude_EE[14648];
  float _amplitude_second_EE[14648];
  float _chi2_EE[14648];
  float _chi2_second_EE[14648];
  float _amplitudeError_EE[14648];
  float _amplitudeError_second_EE[14648];
  float _digi_ped_subtracted_EE[14648*10];
  int   _ix[14648];
  int   _iy[14648];
  int   _iz[14648];

  int   _size_EB;
  int   _size_EE;
  
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
PulseTreeProducer::PulseTreeProducer(const edm::ParameterSet& iConfig) //,  edm::ConsumesCollector& myConsumesCollector)

{
  //now do what ever initialization is needed
  usesResource("TFileService");
  
  
  //now do what ever initialization is needed
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  
  _token_ebrechits = consumes<EcalUncalibratedRecHitCollection>(iConfig.getParameter<edm::InputTag>("EcalUncalibRecHitsEBCollection"));
  _token_eerechits = consumes<EcalUncalibratedRecHitCollection>(iConfig.getParameter<edm::InputTag>("EcalUncalibRecHitsEECollection"));
  
  _token_second_ebrechits = consumes<EcalUncalibratedRecHitCollection>(iConfig.getParameter<edm::InputTag>("EcalUncalibRecHitsEBCollectionSecond"));
  _token_second_eerechits = consumes<EcalUncalibratedRecHitCollection>(iConfig.getParameter<edm::InputTag>("EcalUncalibRecHitsEECollectionSecond"));
  
  _token_ebdigi = consumes<EBDigiCollection>(iConfig.getParameter<edm::InputTag>("EBDigiCollection"));
  _token_eedigi = consumes<EEDigiCollection>(iConfig.getParameter<edm::InputTag>("EEDigiCollection"));
  
  
//   _pedsToken = myConsumesCollector.esConsumes<EcalPedestals, EcalPedestalsRcd>();
  
  _outTree = fs->make<TTree>("tree","tree");
  
  _outTree->Branch("run",               &_run,             "run/i");
  _outTree->Branch("lumi",              &_lumi,            "lumi/s");
  _outTree->Branch("bx",                &_bx,              "bx/s");
  _outTree->Branch("event",             &_event,           "event/i");
  
  _outTree->Branch("digi_ped_subtracted_EB",        _digi_ped_subtracted_EB,        "digi_ped_subtracted_EB[612000]/F"); // 61200*10
  _outTree->Branch("amplitude_EB",        _amplitude_EB,        "amplitude_EB[61200]/F");
  _outTree->Branch("amplitude_second_EB", _amplitude_second_EB, "amplitude_second_EB[61200]/F");
  _outTree->Branch("chi2_EB",             _chi2_EB,             "chi2_EB[61200]/F");
  _outTree->Branch("chi2_second_EB",      _chi2_second_EB,      "chi2_second_EB[61200]/F");
  _outTree->Branch("amplitudeError_EB",             _amplitudeError_EB,             "amplitudeError_EB[61200]/F");
  _outTree->Branch("amplitudeError_second_EB",      _amplitudeError_second_EB,      "amplitudeError_second_EB[61200]/F");
  _outTree->Branch("ieta",                _ieta,                "ieta[61200]/I");
  _outTree->Branch("iphi",                _iphi,                "iphi[61200]/I");
  
  _outTree->Branch("digi_ped_subtracted_EE",        _digi_ped_subtracted_EE,        "digi_ped_subtracted_EE[146480]/F"); // 14648*10
  _outTree->Branch("amplitude_EE",        _amplitude_EE,        "amplitude_EE[14648]/F");
  _outTree->Branch("amplitude_second_EE", _amplitude_second_EE, "amplitude_second_EE[14648]/F");
  _outTree->Branch("chi2_EE",             _chi2_EE,             "chi2_EE[14648]/F");
  _outTree->Branch("chi2_second_EE",      _chi2_second_EE,      "chi2_second_EE[14648]/F");
  _outTree->Branch("amplitudeError_EE",             _amplitudeError_EE,             "amplitudeError_EE[14648]/F");
  _outTree->Branch("amplitudeError_second_EE",      _amplitudeError_second_EE,      "amplitudeError_second_EE[14648]/F");
  _outTree->Branch("ix",                  _ix,                  "ix[14648]/I");
  _outTree->Branch("iy",                  _iy,                  "iy[14648]/I");
  _outTree->Branch("iz",                  _iz,                  "iz[14648]/I");
  
  _outTree->Branch("size_EB",             &_size_EB,                  "size_EB/I");
  _outTree->Branch("size_EE",             &_size_EE,                  "size_EE/I");
  
  
  
}


PulseTreeProducer::~PulseTreeProducer()
{
  
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  
}


//
// member functions
//

// ------------ method called for each event  ------------
void
PulseTreeProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  _run = iEvent.eventAuxiliary().run();
  _lumi = iEvent.eventAuxiliary().luminosityBlock();
  _bx = iEvent.eventAuxiliary().bunchCrossing();
  _event = iEvent.eventAuxiliary().event();
  
  
  
  
  //---- pedestals
  edm::ESHandle< EcalPedestals > ecalPedestals;
  iSetup.get<EcalPedestalsRcd>().get(ecalPedestals);
  _peds = ecalPedestals.product();
  
  
  
  
  
  //---- uncalibrated rechits
  edm::Handle<EcalUncalibratedRecHitCollection> ebrechithandle;
  const EcalUncalibratedRecHitCollection *ebrechits = NULL;
  edm::Handle<EcalUncalibratedRecHitCollection> eerechithandle;
  const EcalUncalibratedRecHitCollection *eerechits = NULL;
  
  iEvent.getByToken(_token_ebrechits,ebrechithandle);
  ebrechits = ebrechithandle.product();
  iEvent.getByToken(_token_eerechits,eerechithandle);
  eerechits = eerechithandle.product();
  
  
  edm::Handle<EcalUncalibratedRecHitCollection> ebrechit_second_handle;
  const EcalUncalibratedRecHitCollection *ebrechits_second = NULL;
  edm::Handle<EcalUncalibratedRecHitCollection> eerechit_second_handle;
  const EcalUncalibratedRecHitCollection *eerechits_second = NULL;
  
  iEvent.getByToken(_token_second_ebrechits,ebrechit_second_handle);
  ebrechits_second = ebrechit_second_handle.product();
  iEvent.getByToken(_token_second_eerechits,eerechit_second_handle);
  eerechits_second = eerechit_second_handle.product();
  
  
  //---- digis
  edm::Handle<EBDigiCollection> ebdigihandle;
  const EBDigiCollection *ebdigis = NULL;
  edm::Handle<EEDigiCollection> eedigihandle;
  const EEDigiCollection *eedigis = NULL;
  
  iEvent.getByToken(_token_ebdigi,ebdigihandle);
  ebdigis = ebdigihandle.product();
  iEvent.getByToken(_token_eedigi,eedigihandle);
  eedigis = eedigihandle.product();
  
  
  
  
  //---- setup default
  for (int ixtal=0; ixtal < 61200; ixtal++) {
    for (int i=0; i<10; i++) _digi_ped_subtracted_EB[ixtal*10+i] = -999;
    _amplitude_EB[ixtal] = -999;
    _amplitude_second_EB[ixtal] = -999;
    _chi2_EB[ixtal] = -999;
    _chi2_second_EB[ixtal] = -999;
    _amplitudeError_EB[ixtal] = -999;
    _amplitudeError_second_EB[ixtal] = -999;
    _ieta[ixtal] = -999;
    _iphi[ixtal] = -999;
  }
  for (int ixtal=0; ixtal < 14648; ixtal++) {
    for (int i=0; i<10; i++) _digi_ped_subtracted_EE[ixtal*10+i] = -999;
    _amplitude_EE[ixtal] = -999;
    _amplitude_second_EE[ixtal] = -999;
    _chi2_EE[ixtal] = -999;
    _chi2_second_EE[ixtal] = -999;
    _amplitudeError_EE[ixtal] = -999;
    _amplitudeError_second_EE[ixtal] = -999;
    _ix[ixtal] = -999;
    _iy[ixtal] = -999;
    _iz[ixtal] = -999;
  }
  
  _size_EB = -999;
  _size_EE = -999;
  
  
  
  //---- geometry 
  
  _size_EB = ebrechits->size();
  //   std::cout << " ebrechits->size() = " << ebrechits->size() << std::endl;
  for (EcalUncalibratedRecHitCollection::const_iterator itrechit = ebrechits->begin(); itrechit != ebrechits->end(); itrechit++ ) {
    _amplitude_EB[EBDetId(itrechit->id()).hashedIndex()] =  itrechit->amplitude();  //----> only in EcalUncalibratedRecHit
    _chi2_EB[EBDetId(itrechit->id()).hashedIndex()] =  itrechit->chi2();  //----> only in EcalUncalibratedRecHit
    _amplitudeError_EB[EBDetId(itrechit->id()).hashedIndex()] =  itrechit->amplitudeError();  //----> only in EcalUncalibratedRecHit
    _ieta[EBDetId(itrechit->id()).hashedIndex()] = EBDetId(itrechit->id()).ieta();
    _iphi[EBDetId(itrechit->id()).hashedIndex()] = EBDetId(itrechit->id()).iphi();    
  }
  
  _size_EE = eerechits->size();
  //   std::cout << " eerechits->size() = " << eerechits->size() << std::endl;
  for (EcalUncalibratedRecHitCollection::const_iterator itrechit = eerechits->begin(); itrechit != eerechits->end(); itrechit++ ) {
    _amplitude_EE[EEDetId(itrechit->id()).hashedIndex()] =  itrechit->amplitude();  //----> only in EcalUncalibratedRecHit
    _chi2_EE[EEDetId(itrechit->id()).hashedIndex()] =  itrechit->chi2();  //----> only in EcalUncalibratedRecHit
    _amplitudeError_EE[EEDetId(itrechit->id()).hashedIndex()] =  itrechit->amplitudeError();  //----> only in EcalUncalibratedRecHit
    _ix[EEDetId(itrechit->id()).hashedIndex()] = EEDetId(itrechit->id()).ix();
    _iy[EEDetId(itrechit->id()).hashedIndex()] = EEDetId(itrechit->id()).iy();
    _iz[EEDetId(itrechit->id()).hashedIndex()] = EEDetId(itrechit->id()).zside();
  }
  
  
  
  for (EcalUncalibratedRecHitCollection::const_iterator itrechit = ebrechits_second->begin(); itrechit != ebrechits_second->end(); itrechit++ ) {
    _amplitude_second_EB[EBDetId(itrechit->id()).hashedIndex()] =  itrechit->amplitude();  //----> only in EcalUncalibratedRecHit
    _chi2_second_EB[EBDetId(itrechit->id()).hashedIndex()] =  itrechit->chi2();  //----> only in EcalUncalibratedRecHit
    _amplitudeError_second_EB[EBDetId(itrechit->id()).hashedIndex()] =  itrechit->amplitudeError();  //----> only in EcalUncalibratedRecHit
    _ieta[EBDetId(itrechit->id()).hashedIndex()] = EBDetId(itrechit->id()).ieta();
    _iphi[EBDetId(itrechit->id()).hashedIndex()] = EBDetId(itrechit->id()).iphi();    
  }
  
  
  for (EcalUncalibratedRecHitCollection::const_iterator itrechit = eerechits_second->begin(); itrechit != eerechits_second->end(); itrechit++ ) {
    _amplitude_second_EE[EEDetId(itrechit->id()).hashedIndex()] =  itrechit->amplitude();  //----> only in EcalUncalibratedRecHit
    _chi2_second_EE[EEDetId(itrechit->id()).hashedIndex()] =  itrechit->chi2();  //----> only in EcalUncalibratedRecHit
    _amplitudeError_second_EE[EEDetId(itrechit->id()).hashedIndex()] =  itrechit->amplitudeError();  //----> only in EcalUncalibratedRecHit
    _ix[EEDetId(itrechit->id()).hashedIndex()] = EEDetId(itrechit->id()).ix();
    _iy[EEDetId(itrechit->id()).hashedIndex()] = EEDetId(itrechit->id()).iy();
    _iz[EEDetId(itrechit->id()).hashedIndex()] = EEDetId(itrechit->id()).zside();
  }
  
  
  
  //---- dump digis
  
  for (EBDigiCollection::const_iterator itdigi = ebdigis->begin(); itdigi != ebdigis->end(); itdigi++ ) {
    
    float pedestal = 0;    
    DetId id = (EBDetId&)((*itdigi));
    pedestal = float((_peds->find(id))->mean_x12);
    
    //                                                           0xFFF = 4095
    for (int iSample = 0; iSample < 10; iSample++) {
      float value = ( int( (*itdigi) [iSample] ) & 0xFFF );
      _digi_ped_subtracted_EB[((EBDetId&)((*itdigi))).hashedIndex() *10 + iSample] = value - pedestal;
    }
  }
  
  
  for (EEDigiCollection::const_iterator itdigi = eedigis->begin(); itdigi != eedigis->end(); itdigi++ ) {
    
    float pedestal = 0;
    DetId id = (EEDetId&)((*itdigi));
    pedestal = float((_peds->find(id))->mean_x12);
    
    //                                                           0xFFF = 4095
    for (int iSample = 0; iSample < 10; iSample++) {
      float value = ( int( (*itdigi) [iSample] ) & 0xFFF );
      _digi_ped_subtracted_EE[((EEDetId&)((*itdigi))).hashedIndex() *10 + iSample] = value - pedestal;
    }
  }
  
  _outTree->Fill();  
  
}




// ------------ method called once each job just before starting event loop  ------------
void 
PulseTreeProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PulseTreeProducer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PulseTreeProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PulseTreeProducer);

