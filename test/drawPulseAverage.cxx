//
// Mean pulse (16 samples, 6.25 ns) over ALL the crystals and ALL the events,
// from the tree of SamplesTreeProducer
//
//   pedestal = (ADC[0] + ADC[1]) / 2, subtracted to all the samples before anything else (also before the cut)
//   (ADC - pedestal) / (ADC[5] - pedestal), only pulses with ADC[5] - pedestal > minAmplitude
//   two canvases saved in outputFolder in .root and .png:
//     pulse_average     mean pulse, error bars = uncertainty on the mean
//     pulse_2d_average  TH2F of all the pulses vs time, with the mean on top
//   and the same with _noNumber, without the number of pulses in the legend
//
// Run (compiled, faster):
//   root -l 'drawPulseAverage.cxx+("samples.root", 50)'
//

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TH2F.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TString.h"
#include "TSystem.h"
#include <iostream>
#include <string>

void drawPulseAverage(std::string inputFile = "samples.root", float minAmplitude = 0, std::string outputFolder = "pulses_output") {

  const int nSamples = 16;
  const float sampleTime = 6.25; // ns
  const int iSampleMax = 5;       // sample at the peak (phase-2 simulation), counting from 0

  gSystem->mkdir(outputFolder.c_str(), true);

  TFile* _file0 = TFile::Open(inputFile.c_str());
  TTree* tree = (TTree*) _file0->Get("SamplesTreeProducer/tree");

  int nDigis;
  static int      ieta[61200];
  static UShort_t adc[61200*16];

  tree->SetBranchStatus("*", 0);
  tree->SetBranchStatus("nDigis", 1);
  tree->SetBranchStatus("nSamplesTot", 1);
  tree->SetBranchStatus("ieta", 1);
  tree->SetBranchStatus("adc", 1);

  tree->SetBranchAddress("nDigis", &nDigis);
  tree->SetBranchAddress("ieta",   ieta);
  tree->SetBranchAddress("adc",    adc);

  TProfile* p_pulse = new TProfile ("p_pulse", "", nSamples, -sampleTime/2., (nSamples-0.5)*sampleTime);
  TH2F*     h_pulse = new TH2F ("h_pulse", "all crystals, all events;time [ns];(ADC - pedestal) / (ADC[5] - pedestal)", nSamples, -sampleTime/2., (nSamples-0.5)*sampleTime, 140, -0.2, 1.2);
  int nPulses = 0;

  //---- loop over events
  for (int iEvent = 0; iEvent < tree->GetEntries(); iEvent++) {
    tree->GetEntry(iEvent);
    if (iEvent % 100 == 0) std::cout << " event " << iEvent << " / " << tree->GetEntries() << std::endl;

    for (int iDigi = 0; iDigi < nDigis; iDigi++) {
      //---- pedestal from the first two samples, subtracted before anything else
      float pedestal = 0.5 * (adc[iDigi*nSamples] + adc[iDigi*nSamples + 1]);
      //---- normalize to the sample at the peak
      float maximum = adc[iDigi*nSamples + iSampleMax] - pedestal;
      if (maximum <= 0 || maximum <= minAmplitude) continue;

      nPulses++;
      for (int iSample = 0; iSample < nSamples; iSample++) {
        float value = (adc[iDigi*nSamples + iSample] - pedestal) / maximum;
        p_pulse->Fill(iSample*sampleTime, value);
        h_pulse->Fill(iSample*sampleTime, value);
      }
    }
  }
  std::cout << " all crystals, all events: " << nPulses << " pulses with ADC[5] - pedestal > " << minAmplitude << " ADC" << std::endl;

  //---- mean and uncertainty on the mean as a graph: red markers and red error bars
  TGraphErrors* gr_mean = new TGraphErrors ();
  std::cout << " mean pulse: " << std::endl << "   ";
  for (int iSample = 0; iSample < nSamples; iSample++) {
    gr_mean->SetPoint(iSample, iSample*sampleTime, p_pulse->GetBinContent(iSample+1));
    gr_mean->SetPointError(iSample, 0, p_pulse->GetBinError(iSample+1));
    std::cout << Form("%.5f", p_pulse->GetBinContent(iSample+1)) << " ";
  }
  std::cout << std::endl;
  gr_mean->SetMarkerStyle(20);
  gr_mean->SetMarkerColor(kRed);
  gr_mean->SetLineColor(kRed);
  gr_mean->SetLineWidth(2);
  gr_mean->SetTitle("all crystals, all events;time [ns];(ADC - pedestal) / (ADC[5] - pedestal)");

  TLegend* leg = new TLegend (0.60, 0.80, 0.88, 0.88);
  leg->SetBorderSize(1);
  leg->AddEntry(gr_mean, Form("mean, %d pulses", nPulses), "lp");

  TCanvas* cc = new TCanvas ("pulse_average", "", 800, 600);
  gr_mean->Draw("AP");
  leg->Draw();
  cc->SaveAs((outputFolder + "/pulse_average.png").c_str());
  cc->SaveAs((outputFolder + "/pulse_average.root").c_str());

  //---- all the pulses (TH2F) with the mean on top
  TCanvas* cc_2d = new TCanvas ("pulse_2d_average", "", 800, 600);
  h_pulse->SetStats(0);
  h_pulse->Draw("colz");
  gr_mean->Draw("P same");
  leg->Draw();
  cc_2d->SaveAs((outputFolder + "/pulse_2d_average.png").c_str());
  cc_2d->SaveAs((outputFolder + "/pulse_2d_average.root").c_str());

  //---- same plots without the number of pulses in the legend
  ((TLegendEntry*) leg->GetListOfPrimitives()->First())->SetLabel("mean");
  //---- smaller box, sized to the text, bigger text, top right
  leg->SetX1NDC(0.72);
  leg->SetX2NDC(0.88);
  leg->SetY1NDC(0.80);
  leg->SetY2NDC(0.88);
  leg->SetTextSize(0.045);
  cc->Modified();
  cc->SaveAs((outputFolder + "/pulse_average_noNumber.png").c_str());
  cc->SaveAs((outputFolder + "/pulse_average_noNumber.root").c_str());
  cc_2d->Modified();
  cc_2d->SaveAs((outputFolder + "/pulse_2d_average_noNumber.png").c_str());
  cc_2d->SaveAs((outputFolder + "/pulse_2d_average_noNumber.root").c_str());

}
