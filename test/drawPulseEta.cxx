//
// Draw the mean pulse (16 samples, 6.25 ns) of one ieta ring, all the crystals along iphi,
// from the tree of SamplesTreeProducer
//
//   all events and all iphi overlaid, ADC / ADC[5], only pulses with ADC[5] > minAmplitude
//   red points = mean pulse, error bars = uncertainty on the mean, also printed per sample
//
// Run:
//   root -l 'drawPulseEta.cxx("samples.root", 10, 50)'
//

void drawPulseEta(std::string inputFile = "samples.root", int ieta_sel = 1, float minAmplitude = 0) {

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

  TString title = Form("ieta = %d, all iphi", ieta_sel);

  TCanvas* cc = new TCanvas ("cc", "", 800, 600);

  //---- all events and all iphi overlaid, normalized to the sample at the peak
  TProfile* p_pulse = new TProfile ("p_pulse", "", nSamples, -sampleTime/2., (nSamples-0.5)*sampleTime);

  int nPulses = 0;
  std::set<int> iphi_used;
  for (int iEvent = 0; iEvent < tree->GetEntries(); iEvent++) {
    tree->GetEntry(iEvent);

    for (int iDigi = 0; iDigi < nDigis; iDigi++) {
      if (ieta[iDigi] != ieta_sel) continue;

      //---- normalize to the sample at the peak, no pedestal subtraction
      float maximum = adc[iDigi*nSamples + iSampleMax];
      if (maximum <= 0 || maximum <= minAmplitude) continue;

      nPulses++;
      iphi_used.insert(iphi[iDigi]);
      for (int iSample = 0; iSample < nSamples; iSample++) {
        float value = adc[iDigi*nSamples + iSample] / maximum;
        p_pulse->Fill(iSample*sampleTime, value);
      }
    }
  }
  std::cout << " " << title << ": " << nPulses << " pulses with ADC[5] > " << minAmplitude << " ADC, from " << iphi_used.size() << " crystals" << std::endl;

  //---- mean and uncertainty on the mean per sample
  std::cout << " sample   time [ns]   mean     error" << std::endl;
  for (int iSample = 0; iSample < nSamples; iSample++) {
    std::cout << Form(" %4d   %8.2f   %7.4f  %7.4f", iSample, iSample*sampleTime, p_pulse->GetBinContent(iSample+1), p_pulse->GetBinError(iSample+1)) << std::endl;
  }

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

  TLegend* leg = new TLegend (0.60, 0.80, 0.88, 0.88);
  leg->AddEntry(gr_mean, Form("mean, %d pulses", nPulses), "lp");
  leg->Draw();

}
