//
// Draw the pulse (16 samples, 6.25 ns) of one crystal from the tree of SamplesTreeProducer
//
//   iEvent >= 0: raw ADC of that event (full markers = gain x10, open markers = gain x1)
//   iEvent  < 0: all events overlaid, ADC / ADC[5], only pulses with ADC[5] > minAmplitude
//                red points = mean pulse, error bars = uncertainty on the mean, also printed per sample
//                second canvas: TH2F of all the pulses, ADC / ADC[5] vs time, with the mean on top
//
// Run:
//   root -l 'drawPulse.cxx("samples.root", 10, 100, 0)'
//   root -l 'drawPulse.cxx("samples.root", 10, 100, -1, 50)'
//

void drawPulse(std::string inputFile = "samples.root", int ieta_sel = 1, int iphi_sel = 1, int iEvent_sel = -1, float minAmplitude = 0) {

  const int nSamples = 16;
  const float sampleTime = 6.25; // ns
  const int iSampleMax = 5;       // sample at the peak (phase-2 simulation), counting from 0

  TFile* _file0 = TFile::Open(inputFile.c_str());
  TTree* tree = (TTree*) _file0->Get("SamplesTreeProducer/tree");

  int nDigis;
  static int      ieta[61200];
  static int      iphi[61200];
  static UShort_t adc[61200*16];
  static UChar_t  gain[61200*16];

  tree->SetBranchStatus("*", 0);
  tree->SetBranchStatus("nDigis", 1);
  tree->SetBranchStatus("nSamplesTot", 1);
  tree->SetBranchStatus("ieta", 1);
  tree->SetBranchStatus("iphi", 1);
  tree->SetBranchStatus("adc", 1);
  tree->SetBranchStatus("gain", 1);

  tree->SetBranchAddress("nDigis", &nDigis);
  tree->SetBranchAddress("ieta",   ieta);
  tree->SetBranchAddress("iphi",   iphi);
  tree->SetBranchAddress("adc",    adc);
  tree->SetBranchAddress("gain",   gain);

  TString title = Form("ieta = %d, iphi = %d", ieta_sel, iphi_sel);

  TCanvas* cc = new TCanvas ("cc", "", 800, 600);

  //---- single event
  if (iEvent_sel >= 0) {
    tree->GetEntry(iEvent_sel);

    TGraph* gr_g10 = new TGraph ();
    TGraph* gr_g1  = new TGraph ();
    bool found = false;
    for (int iDigi = 0; iDigi < nDigis; iDigi++) {
      if (ieta[iDigi] != ieta_sel || iphi[iDigi] != iphi_sel) continue;
      found = true;
      for (int iSample = 0; iSample < nSamples; iSample++) {
        int index = iDigi*nSamples + iSample;
        if (gain[index] == 0) gr_g10->SetPoint(gr_g10->GetN(), iSample*sampleTime, adc[index]);
        else                  gr_g1 ->SetPoint(gr_g1 ->GetN(), iSample*sampleTime, adc[index]);
      }
    }
    if (!found) {
      std::cout << " crystal " << title << " not found in event " << iEvent_sel << std::endl;
      return;
    }

    TMultiGraph* mg = new TMultiGraph ();
    gr_g10->SetMarkerStyle(20);
    gr_g10->SetMarkerColor(kRed);
    gr_g1 ->SetMarkerStyle(24);
    gr_g1 ->SetMarkerColor(kRed);
    if (gr_g10->GetN()) mg->Add(gr_g10, "P");
    if (gr_g1 ->GetN()) mg->Add(gr_g1,  "P");
    mg->SetTitle(title + Form(", event %d;time [ns];ADC", iEvent_sel));
    mg->Draw("A");
    return;
  }

  //---- all events overlaid, normalized to the sample at the peak
//   TH2F* h_pulse = new TH2F ("h_pulse", title + ";time [ns];ADC / ADC[5]", nSamples, -sampleTime/2., (nSamples-0.5)*sampleTime, 120, -0.1, 1.1);
//   TProfile* p_pulse = new TProfile ("p_pulse", "", nSamples, -sampleTime/2., (nSamples-0.5)*sampleTime, "s");
  TProfile* p_pulse = new TProfile ("p_pulse", "", nSamples, -sampleTime/2., (nSamples-0.5)*sampleTime);
  TH2F* h_pulse = new TH2F ("h_pulse", title + ";time [ns];ADC / ADC[5]", nSamples, -sampleTime/2., (nSamples-0.5)*sampleTime, 140, -0.2, 1.2);

  int nPulses = 0;
  for (int iEvent = 0; iEvent < tree->GetEntries(); iEvent++) {
    tree->GetEntry(iEvent);

    for (int iDigi = 0; iDigi < nDigis; iDigi++) {
      if (ieta[iDigi] != ieta_sel || iphi[iDigi] != iphi_sel) continue;

      float pedestal = adc[iDigi*nSamples];
      //---- normalize to the sample at the peak, no pedestal subtraction
      float maximum = adc[iDigi*nSamples + iSampleMax];
//       float maximum = adc[iDigi*nSamples + iSampleMax] - pedestal;
//       float maximum = 0;
//       for (int iSample = 0; iSample < nSamples; iSample++) {
//         float value = adc[iDigi*nSamples + iSample] - pedestal;
//         if (value > maximum) maximum = value;
//       }
      if (maximum <= 0 || maximum <= minAmplitude) continue;
//       if (maximum <= 0 || maximum - pedestal <= minAmplitude) continue;

      nPulses++;
      for (int iSample = 0; iSample < nSamples; iSample++) {
        float value = adc[iDigi*nSamples + iSample] / maximum;
//         float value = (adc[iDigi*nSamples + iSample] - pedestal) / maximum;
//         h_pulse->Fill(iSample*sampleTime, value);
        p_pulse->Fill(iSample*sampleTime, value);
        h_pulse->Fill(iSample*sampleTime, value);
      }
    }
  }
  std::cout << " " << title << ": " << nPulses << " pulses with ADC[5] > " << minAmplitude << " ADC" << std::endl;

  //---- mean and uncertainty on the mean per sample
  std::cout << " sample   time [ns]   mean     error" << std::endl;
  for (int iSample = 0; iSample < nSamples; iSample++) {
    std::cout << Form(" %4d   %8.2f   %7.4f  %7.4f", iSample, iSample*sampleTime, p_pulse->GetBinContent(iSample+1), p_pulse->GetBinError(iSample+1)) << std::endl;
  }

//   h_pulse->SetStats(0);
//   h_pulse->Draw("colz");

  //---- mean and uncertainty on the mean as a graph: red markers and red error bars
  TGraphErrors* gr_mean = new TGraphErrors ();
  for (int iSample = 0; iSample < nSamples; iSample++) {
    gr_mean->SetPoint(iSample, iSample*sampleTime, p_pulse->GetBinContent(iSample+1));
    gr_mean->SetPointError(iSample, 0, p_pulse->GetBinError(iSample+1));
  }
  gr_mean->SetMarkerStyle(20);
  gr_mean->SetMarkerColor(kRed);
  gr_mean->SetLineColor(kRed);
  gr_mean->SetLineWidth(2);
  gr_mean->SetTitle(title + ";time [ns];ADC / ADC[5]");
  gr_mean->Draw("AP");
//   gr_mean->Draw("P same");

  TLegend* leg = new TLegend (0.60, 0.80, 0.88, 0.88);
  leg->AddEntry(gr_mean, Form("mean, %d pulses", nPulses), "lp");
  leg->Draw();

  //---- second canvas: all the pulses (TH2F) with the mean on top
  TCanvas* cc_2d = new TCanvas ("cc_2d", "", 800, 600);
  h_pulse->SetStats(0);
  h_pulse->Draw("colz");
  gr_mean->Draw("P same");
  leg->Draw();

}
