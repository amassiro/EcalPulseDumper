//
// Mean pulse (16 samples, 6.25 ns) for every ieta ring, all the crystals along iphi,
// from the tree of SamplesTreeProducer, in one single loop over the tree
//
//   pedestal = (ADC[0] + ADC[1]) / 2, subtracted to all the samples before anything else (also before the cut)
//   all events and all iphi overlaid, (ADC - pedestal) / (ADC[5] - pedestal), only pulses with ADC[5] - pedestal > minAmplitude
//   for each ieta (-85 ... 85, no 0) two canvases saved in outputFolder in .root and .png:
//     pulse_ieta_X     mean pulse, error bars = uncertainty on the mean
//     pulse_2d_ieta_X  TH2F of all the pulses vs time, with the mean on top
//   (X = m85 ... m1, 1 ... 85)
//
// Run (compiled, faster):
//   root -l -b -q 'drawPulseAllEta.cxx+("samples.root", 50)'
// or interpreted:
//   root -l -b -q 'drawPulseAllEta.cxx("samples.root", 50)'
//

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TH2F.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include <iostream>
#include <string>

void drawPulseAllEta(std::string inputFile = "samples.root", float minAmplitude = 0, std::string outputFolder = "pulses_output") {

  const int nSamples = 16;
  const float sampleTime = 6.25; // ns
  const int iSampleMax = 5;       // sample at the peak (phase-2 simulation), counting from 0
  const int nEta = 85;            // ieta = -85 ... 85, no 0

  gROOT->SetBatch(true);
  gSystem->mkdir(outputFolder.c_str(), true);

  TFile* _file0 = TFile::Open(inputFile.c_str());
  TTree* tree = (TTree*) _file0->Get("SamplesTreeProducer/tree");

  int nDigis;
  static int      ieta[61200];
  static int      iphi[61200];
  static UShort_t adc[61200*16];

  tree->SetBranchStatus("*", 0);
  tree->SetBranchStatus("nDigis", 1);
  tree->SetBranchStatus("nSamplesTot", 1);
  tree->SetBranchStatus("ieta", 1);
  tree->SetBranchStatus("iphi", 1);
  tree->SetBranchStatus("adc", 1);

  tree->SetBranchAddress("nDigis", &nDigis);
  tree->SetBranchAddress("ieta",   ieta);
  tree->SetBranchAddress("iphi",   iphi);
  tree->SetBranchAddress("adc",    adc);

  //---- one TProfile and one TH2F per ieta, index = ieta + nEta
  TProfile* p_pulse[2*nEta+1];
  TH2F*     h_pulse[2*nEta+1];
  int       nPulses[2*nEta+1];
  for (int iEta = -nEta; iEta <= nEta; iEta++) {
    int index = iEta + nEta;
    p_pulse[index] = nullptr;
    h_pulse[index] = nullptr;
    nPulses[index] = 0;
    if (iEta == 0) continue;
    TString title = Form("ieta = %d, all iphi", iEta);
    p_pulse[index] = new TProfile (Form("p_pulse_%d", index), "", nSamples, -sampleTime/2., (nSamples-0.5)*sampleTime);
    h_pulse[index] = new TH2F (Form("h_pulse_%d", index), title + ";time [ns];(ADC - pedestal) / (ADC[5] - pedestal)", nSamples, -sampleTime/2., (nSamples-0.5)*sampleTime, 140, -0.2, 1.2);
  }

  //---- loop over events, once
  for (int iEvent = 0; iEvent < tree->GetEntries(); iEvent++) {
    tree->GetEntry(iEvent);
    if (iEvent % 100 == 0) std::cout << " event " << iEvent << " / " << tree->GetEntries() << std::endl;

    for (int iDigi = 0; iDigi < nDigis; iDigi++) {
      int index = ieta[iDigi] + nEta;
      if (index < 0 || index > 2*nEta || !p_pulse[index]) continue;

      //---- pedestal from the first two samples, subtracted before anything else
      float pedestal = 0.5 * (adc[iDigi*nSamples] + adc[iDigi*nSamples + 1]);
      //---- normalize to the sample at the peak
      float maximum = adc[iDigi*nSamples + iSampleMax] - pedestal;
      if (maximum <= 0 || maximum <= minAmplitude) continue;

      nPulses[index]++;
      for (int iSample = 0; iSample < nSamples; iSample++) {
        float value = (adc[iDigi*nSamples + iSample] - pedestal) / maximum;
        p_pulse[index]->Fill(iSample*sampleTime, value);
        h_pulse[index]->Fill(iSample*sampleTime, value);
      }
    }
  }

  //---- draw and save
  for (int iEta = -nEta; iEta <= nEta; iEta++) {
    if (iEta == 0) continue;
    int index = iEta + nEta;
    TString title = Form("ieta = %d, all iphi", iEta);
    TString name  = (iEta < 0) ? Form("ieta_m%d", -iEta) : Form("ieta_%d", iEta);
    std::cout << " " << title << ": " << nPulses[index] << " pulses with ADC[5] - pedestal > " << minAmplitude << " ADC" << std::endl;

    //---- mean and uncertainty on the mean as a graph: red markers and red error bars
    TGraphErrors* gr_mean = new TGraphErrors ();
    for (int iSample = 0; iSample < nSamples; iSample++) {
      gr_mean->SetPoint(iSample, iSample*sampleTime, p_pulse[index]->GetBinContent(iSample+1));
      gr_mean->SetPointError(iSample, 0, p_pulse[index]->GetBinError(iSample+1));
    }
    gr_mean->SetMarkerStyle(20);
    gr_mean->SetMarkerColor(kRed);
    gr_mean->SetLineColor(kRed);
    gr_mean->SetLineWidth(2);
    gr_mean->SetTitle(title + ";time [ns];(ADC - pedestal) / (ADC[5] - pedestal)");

    TLegend* leg = new TLegend (0.60, 0.80, 0.88, 0.88);
    leg->AddEntry(gr_mean, Form("mean, %d pulses", nPulses[index]), "lp");

    TCanvas* cc = new TCanvas ("pulse_" + name, "", 800, 600);
    gr_mean->Draw("AP");
    leg->Draw();
    cc->SaveAs((outputFolder + "/pulse_" + name + ".png").Data());
    cc->SaveAs((outputFolder + "/pulse_" + name + ".root").Data());

    //---- all the pulses (TH2F) with the mean on top
    TCanvas* cc_2d = new TCanvas ("pulse_2d_" + name, "", 800, 600);
    h_pulse[index]->SetStats(0);
    h_pulse[index]->Draw("colz");
    gr_mean->Draw("P same");
    leg->Draw();
    cc_2d->SaveAs((outputFolder + "/pulse_2d_" + name + ".png").Data());
    cc_2d->SaveAs((outputFolder + "/pulse_2d_" + name + ".root").Data());

    delete cc;
    delete cc_2d;
  }

}
