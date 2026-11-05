//
// Apply the amplitude weights per ieta (from LocalCreateWeightsFromPulses.cxx) to one pulse
// of one crystal (ieta, iphi) from the tree of SamplesTreeProducer
//
//   weights file: "# ieta weights", then one line per ieta: ieta w[0] ... w[15]
//   pedestal  = (ADC[0] + ADC[1]) / 2
//   amplitude = sum_i w[i] * ADC[i]
//   event: iEvent >= 0 that event, iEvent < 0 first event with ADC[5] - pedestal > threshold
//   then buttons "-" and "+" in the canvas to show the previous / next event with ADC[5] - pedestal > threshold
//
// Plot: ADC vs time, dashed line = pedestal, blue square = amplitude + pedestal (next to sample 5)
//
// Run:
//   root -l 'ApplyWeights.cxx("samples.root", "weights_pulses_all.txt", -54, 119, 50)'
//   root -l 'ApplyWeights.cxx("samples.root", "weights_pulses_all.txt", -54, 119, 50, 3)'
//

//---- global, used by the buttons
const int aw_nSamples = 16;
const float aw_sampleTime = 6.25; // ns
const int aw_iSampleMax = 5;       // sample at the peak (phase-2 simulation), counting from 0

TTree* aw_tree = nullptr;
TCanvas* aw_cc = nullptr;
TPad* aw_pad = nullptr;
std::vector<float> aw_weights;
int aw_ieta_sel;
int aw_iphi_sel;
float aw_threshold;
int aw_iEventCurrent = -1;

int aw_nDigis;
int      aw_ieta[61200];
int      aw_iphi[61200];
UShort_t aw_adc[61200*16];
UChar_t  aw_gain[61200*16];


//---- find the digi of the selected crystal with a pulse, from iEventStart forward (direction = 1) or backward (direction = -1)
//---- onlyOneEvent = look only at iEventStart. Return the event, -1 if not found
int ApplyWeightsFind(int iEventStart, int direction, bool onlyOneEvent, int& iDigiFound) {
  int iEventEnd = onlyOneEvent ? iEventStart + direction : (direction > 0 ? aw_tree->GetEntries() : -1);
  for (int iEvent = iEventStart; iEvent != iEventEnd && iEvent >= 0 && iEvent < aw_tree->GetEntries(); iEvent += direction) {
    aw_tree->GetEntry(iEvent);
    for (int iDigi = 0; iDigi < aw_nDigis; iDigi++) {
      if (aw_ieta[iDigi] != aw_ieta_sel || aw_iphi[iDigi] != aw_iphi_sel) continue;
      float pedestal = 0.5 * (aw_adc[iDigi*aw_nSamples] + aw_adc[iDigi*aw_nSamples + 1]);
      if (aw_adc[iDigi*aw_nSamples + aw_iSampleMax] - pedestal > aw_threshold) {
        iDigiFound = iDigi;
        return iEvent;
      }
      break;
    }
  }
  return -1;
}


//---- draw the pulse of event iEvent (already found with ApplyWeightsFind, tree entry loaded)
void ApplyWeightsDraw(int iEvent, int iDigi) {

  const int nSamples = aw_nSamples;
  TString title = Form("ieta = %d, iphi = %d", aw_ieta_sel, aw_iphi_sel);

  //---- pedestal and amplitude with the weights
  float pedestal = 0.5 * (aw_adc[iDigi*nSamples] + aw_adc[iDigi*nSamples + 1]);
  float amplitude = 0;
  for (int iSample = 0; iSample < nSamples; iSample++) amplitude += aw_weights[iSample] * aw_adc[iDigi*nSamples + iSample];

  std::cout << " " << title << ", event " << iEvent << std::endl;
  std::cout << "   pedestal           = " << pedestal << std::endl;
  std::cout << "   ADC[5] - pedestal  = " << aw_adc[iDigi*nSamples + aw_iSampleMax] - pedestal << std::endl;
  std::cout << "   amplitude weights  = " << amplitude << std::endl;
  std::cout << "   gain               = " << (aw_gain[iDigi*nSamples + aw_iSampleMax] == 0 ? "x10" : "x1") << std::endl;

  //---- pulse
  TGraph* gr_pulse = new TGraph ();
  for (int iSample = 0; iSample < nSamples; iSample++) gr_pulse->SetPoint(iSample, iSample*aw_sampleTime, aw_adc[iDigi*nSamples + iSample]);
  gr_pulse->SetMarkerStyle(20);
  gr_pulse->SetMarkerColor(kRed);
  gr_pulse->SetTitle(title + Form(", event %d;time [ns];ADC", iEvent));
  gr_pulse->SetMinimum(0);

  //---- amplitude + pedestal
  TGraph* gr_amplitude = new TGraph ();
  //---- slightly shifted to the right of sample 5, to be visible
  gr_amplitude->SetPoint(0, aw_iSampleMax*aw_sampleTime + 3., amplitude + pedestal);
  gr_amplitude->SetMarkerStyle(21);
  gr_amplitude->SetMarkerSize(1.5);
  gr_amplitude->SetMarkerColor(kBlue);

  aw_pad->Clear();
  aw_pad->cd();
  gr_pulse->Draw("AP");
  gr_amplitude->Draw("P same");

  //---- pedestal, dashed line
  TLine* line_pedestal = new TLine (gr_pulse->GetXaxis()->GetXmin(), pedestal, gr_pulse->GetXaxis()->GetXmax(), pedestal);
  line_pedestal->SetLineStyle(2);
  line_pedestal->SetLineWidth(2);
  line_pedestal->Draw();

  TLegend* leg = new TLegend (0.55, 0.75, 0.88, 0.88);
  leg->AddEntry(gr_pulse, "ADC", "p");
  leg->AddEntry(gr_amplitude, Form("amplitude + pedestal = %.1f", amplitude + pedestal), "p");
  leg->AddEntry(line_pedestal, Form("pedestal = %.1f", pedestal), "l");
  leg->Draw();

  aw_pad->Modified();
  aw_cc->Modified();
  aw_cc->Update();
}


//---- called by the buttons: direction = 1 next event, -1 previous event
void ApplyWeightsMove(int direction) {
  int iDigi = -1;
  int iEvent = ApplyWeightsFind(aw_iEventCurrent + direction, direction, false, iDigi);
  if (iEvent < 0) {
    std::cout << " no " << (direction > 0 ? "next" : "previous") << " pulse with ADC[5] - pedestal > " << aw_threshold << std::endl;
    return;
  }
  aw_iEventCurrent = iEvent;
  ApplyWeightsDraw(iEvent, iDigi);
}


void ApplyWeights(std::string inputFile = "samples.root", std::string weightsFile = "weights_pulses_all.txt", int ieta_sel = 1, int iphi_sel = 1, float threshold = 50, int iEvent_sel = -1) {

  aw_ieta_sel = ieta_sel;
  aw_iphi_sel = iphi_sel;
  aw_threshold = threshold;

  //---- read the weights of the selected ieta
  aw_weights.clear();
  std::ifstream file (weightsFile.c_str());
  if (!file.is_open()) {
    std::cout << " Error: cannot open " << weightsFile << std::endl;
    return;
  }
  std::string line;
  while (std::getline(file, line)) {
    if (line.size() == 0 || line[0] == '#') continue;
    std::istringstream iss (line);
    int ieta;
    iss >> ieta;
    if (ieta != ieta_sel) continue;
    float value;
    while (iss >> value) aw_weights.push_back(value);
  }
  file.close();
  if (aw_weights.size() != aw_nSamples) {
    std::cout << " Error: weights for ieta = " << ieta_sel << " not found in " << weightsFile << std::endl;
    return;
  }

  TFile* _file0 = TFile::Open(inputFile.c_str());
  aw_tree = (TTree*) _file0->Get("SamplesTreeProducer/tree");

  aw_tree->SetBranchStatus("*", 0);
  aw_tree->SetBranchStatus("nDigis", 1);
  aw_tree->SetBranchStatus("nSamplesTot", 1);
  aw_tree->SetBranchStatus("ieta", 1);
  aw_tree->SetBranchStatus("iphi", 1);
  aw_tree->SetBranchStatus("adc", 1);
  aw_tree->SetBranchStatus("gain", 1);

  aw_tree->SetBranchAddress("nDigis", &aw_nDigis);
  aw_tree->SetBranchAddress("ieta",   aw_ieta);
  aw_tree->SetBranchAddress("iphi",   aw_iphi);
  aw_tree->SetBranchAddress("adc",    aw_adc);
  aw_tree->SetBranchAddress("gain",   aw_gain);

  //---- first event
  int iDigi = -1;
  int iEvent = ApplyWeightsFind((iEvent_sel >= 0) ? iEvent_sel : 0, 1, (iEvent_sel >= 0), iDigi);
  if (iEvent < 0) {
    std::cout << " no pulse with ADC[5] - pedestal > " << threshold << " for ieta = " << ieta_sel << ", iphi = " << iphi_sel << std::endl;
    return;
  }
  aw_iEventCurrent = iEvent;

  //---- canvas: pad for the plot, buttons below
  aw_cc = new TCanvas ("cc", "", 800, 650);
  aw_pad = new TPad ("pad", "", 0., 0.08, 1., 1.);
  aw_pad->Draw();

  aw_cc->cd();
  TButton* button_previous = new TButton ("-", "ApplyWeightsMove(-1)", 0.40, 0.01, 0.49, 0.07);
  button_previous->SetTextSize(0.6);
  button_previous->Draw();
  TButton* button_next = new TButton ("+", "ApplyWeightsMove(1)", 0.51, 0.01, 0.60, 0.07);
  button_next->SetTextSize(0.6);
  button_next->Draw();

  ApplyWeightsDraw(iEvent, iDigi);

}


//
// Apply phase-2 weights (from EcalWeights/LocalCreateWeightsPhase2.cxx) to the tree of SamplesTreeProducer
//
//   amplitude = sum_i w_amp[i] * adc[i] * gain_ratio[i]
//   time [ns] = (sum_i w_time[i] * adc[i] * gain_ratio[i]) / amplitude    (dt > 0 = pulse earlier)
//
// gain_ratio = gains[0] / gains[gainId], gains = {10, 1} (ecalPh2::gains): amplitude in gain x10 ADC counts
// (the simulation LSB is Emax / (4095 * gains[gainId]), see EcalLiteDTUCoder.cc).
// In phase-2, when gain x10 saturates the whole frame is read in gain x1.
//
// If amplitude_cmssw is in the tree, it is compared with the local amplitude computed with the
// EcalUncalibRecHitPhase2WeightsProducer convention: adc * gains[gainId]
//
// NB: old version, weights file in the LocalCreateWeightsPhase2.cxx format (one row of 16 weights per parameter)
//
// Run:
//   root -l
//   .L ApplyWeights.cxx
//   ApplyWeightsHistos("samples.root", "weights_TB2025_v1.txt")
//

void ApplyWeightsHistos(std::string inputFile = "samples.root", std::string weightsFile = "weights_Phase2.txt", float minAmplitude = 50) {

  const int nSamples = 16;
  float gains[2] = {10., 1.};

  //---- read weights: first row amplitude, second row amplitude*dt (if in the "# rows:" header)
  std::vector < std::vector<float> > weights;
  bool doTime = false;
  std::ifstream file (weightsFile.c_str());
  std::string line;
  while (std::getline(file, line)) {
    if (line.find("# rows:") == 0 && line.find("amplitude*dt") != std::string::npos) doTime = true;
    if (line.size() == 0 || line[0] == '#') continue;
    std::istringstream iss (line);
    std::vector<float> row;
    float value;
    while (iss >> value) row.push_back(value);
    if (row.size() == nSamples) weights.push_back(row);
  }
  file.close();
  std::cout << " read " << weights.size() << " rows of weights from " << weightsFile << std::endl;
  if (weights.size() == 0) return;
  if (doTime && weights.size() < 2) doTime = false;
  std::cout << " time weights: " << doTime << std::endl;

  TFile* _file0 = TFile::Open(inputFile.c_str());
  TTree* tree = (TTree*) _file0->Get("SamplesTreeProducer/tree");

  int nDigis;
  int nSamplesTot;
  static UShort_t adc[61200*16];
  static UChar_t  gain[61200*16];
  static float    amplitude_cmssw[61200];
  bool hasCMSSW = (tree->GetBranch("amplitude_cmssw") != nullptr);

  tree->SetBranchAddress("nDigis",      &nDigis);
  tree->SetBranchAddress("nSamplesTot", &nSamplesTot);
  tree->SetBranchAddress("adc",         adc);
  tree->SetBranchAddress("gain",        gain);
  if (hasCMSSW) tree->SetBranchAddress("amplitude_cmssw", amplitude_cmssw);

  TH1F* h_amplitude = new TH1F ("h_amplitude", "amplitude [ADC gain x10]", 200, 0, 20000);
  TH1F* h_time      = new TH1F ("h_time", "time [ns]", 200, -10, 10);
  TH1F* h_ratio     = new TH1F ("h_ratio", "amplitude local / CMSSW", 200, 0.8, 1.2);

  //---- loop over events
  for (int iEvent = 0; iEvent < tree->GetEntries(); iEvent++) {
    tree->GetEntry(iEvent);

    for (int iDigi = 0; iDigi < nDigis; iDigi++) {
      float amplitude = 0;
      float amplitude_dt = 0;
      float amplitude_cmssw_convention = 0;
      for (int iSample = 0; iSample < nSamples; iSample++) {
        int gainId = gain[iDigi*nSamples + iSample];
        float value = adc[iDigi*nSamples + iSample] * gains[0] / gains[gainId];
        amplitude += weights[0][iSample] * value;
        if (doTime) amplitude_dt += weights[1][iSample] * value;
        amplitude_cmssw_convention += weights[0][iSample] * adc[iDigi*nSamples + iSample] * gains[gainId];
      }

      h_amplitude->Fill(amplitude);
      if (amplitude < minAmplitude) continue;
      if (doTime) h_time->Fill(amplitude_dt / amplitude);
      if (hasCMSSW && amplitude_cmssw[iDigi] > 0) h_ratio->Fill(amplitude_cmssw_convention / amplitude_cmssw[iDigi]);
    }
  }

  TCanvas* cc = new TCanvas ("cc", "", 1200, 400);
  cc->Divide(3,1);
  cc->cd(1); h_amplitude->Draw(); gPad->SetLogy();
  cc->cd(2); h_time->Draw();
  cc->cd(3); h_ratio->Draw();

}
