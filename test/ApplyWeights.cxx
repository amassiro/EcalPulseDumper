//
// Apply the amplitude weights per ieta (from LocalCreateWeightsFromPulses.cxx) to the pulses
// from the tree of SamplesTreeProducer, one pulse at a time
//
//   weights file: "# ieta weights", then one line per ieta: ieta w[0] ... w[15]
//   pedestal  = (ADC[0] + ADC[1]) / 2
//   amplitude = sum_i w[i] * ADC[i]
//   pulses: ADC[5] - pedestal > threshold
//
//   crystal selection, 0 = not fixed (sequential):
//     ieta = 0, iphi = 0  -->  all crystals, ieta from -85 to 85, and for each ieta iphi from 1 to 360
//     ieta fixed          -->  that ieta, iphi sequential
//     ieta = 0, iphi fixed -->  that iphi, ieta sequential
//   for each crystal all its events with a pulse, then the next crystal
//   iEvent >= 0: start from the first pulse in that event
//
//   buttons "-" and "+" in a separate control bar (TControlBar, as in the TMVA GUI)
//   to show the previous / next pulse
//
// Plot: ADC vs time, dashed line = pedestal, blue square = amplitude + pedestal (next to sample 5)
//
// NB: all the pulses are found in one loop over the tree at the beginning,
//     keep the threshold high enough (with threshold ~ 0 all the 61200 crystals x events are kept)
//
// Run:
//   root -l 'ApplyWeights.cxx("samples.root", "weights_pulses_all.txt", 50)'
//   root -l 'ApplyWeights.cxx("samples.root", "weights_pulses_all.txt", 50, -54)'
//   root -l 'ApplyWeights.cxx("samples.root", "weights_pulses_all.txt", 50, 0, 119)'
//   root -l 'ApplyWeights.cxx("samples.root", "weights_pulses_all.txt", 50, -54, 119, 3)'
//

//---- global, used by the buttons
const int aw_nSamples = 16;
const float aw_sampleTime = 6.25; // ns
const int aw_iSampleMax = 5;       // sample at the peak (phase-2 simulation), counting from 0

TTree* aw_tree = nullptr;
TCanvas* aw_cc = nullptr;
std::map<int, std::vector<float> > aw_weights;
float aw_threshold;

//---- list of the pulses, ordered by ieta, iphi, event
struct aw_pulse {
  int ieta;
  int iphi;
  int iEvent;
  int iDigi;
};
std::vector<aw_pulse> aw_pulses;
int aw_iPulseCurrent = -1;

int aw_nDigis;
int      aw_ieta[61200];
int      aw_iphi[61200];
UShort_t aw_adc[61200*16];
UChar_t  aw_gain[61200*16];


//---- draw the pulse number iPulse of the list
void ApplyWeightsDraw(int iPulse) {

  const int nSamples = aw_nSamples;
  aw_pulse pulse = aw_pulses[iPulse];
  aw_tree->GetEntry(pulse.iEvent);
  int iDigi = pulse.iDigi;
  TString title = Form("ieta = %d, iphi = %d", pulse.ieta, pulse.iphi);

  //---- pedestal and amplitude with the weights of this ieta
  std::vector<float>& weights = aw_weights[pulse.ieta];
  float pedestal = 0.5 * (aw_adc[iDigi*nSamples] + aw_adc[iDigi*nSamples + 1]);
  float amplitude = 0;
  for (int iSample = 0; iSample < nSamples; iSample++) amplitude += weights[iSample] * aw_adc[iDigi*nSamples + iSample];

  std::cout << " pulse " << iPulse + 1 << " / " << aw_pulses.size() << ": " << title << ", event " << pulse.iEvent << std::endl;
  std::cout << "   pedestal           = " << pedestal << std::endl;
  std::cout << "   ADC[5] - pedestal  = " << aw_adc[iDigi*nSamples + aw_iSampleMax] - pedestal << std::endl;
  std::cout << "   amplitude weights  = " << amplitude << std::endl;
  std::cout << "   gain               = " << (aw_gain[iDigi*nSamples + aw_iSampleMax] == 0 ? "x10" : "x1") << std::endl;

  //---- pulse
  TGraph* gr_pulse = new TGraph ();
  for (int iSample = 0; iSample < nSamples; iSample++) gr_pulse->SetPoint(iSample, iSample*aw_sampleTime, aw_adc[iDigi*nSamples + iSample]);
  gr_pulse->SetMarkerStyle(20);
  gr_pulse->SetMarkerColor(kRed);
  gr_pulse->SetTitle(title + Form(", event %d;time [ns];ADC", pulse.iEvent));
  gr_pulse->SetMinimum(0);

  //---- amplitude + pedestal
  TGraph* gr_amplitude = new TGraph ();
  //---- slightly shifted to the right of sample 5, to be visible
  gr_amplitude->SetPoint(0, aw_iSampleMax*aw_sampleTime + 3., amplitude + pedestal);
  gr_amplitude->SetMarkerStyle(21);
  gr_amplitude->SetMarkerSize(1.5);
  gr_amplitude->SetMarkerColor(kBlue);

  aw_cc->Clear();
  aw_cc->cd();
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

  aw_cc->Modified();
  aw_cc->Update();
}


//---- called by the buttons: direction = 1 next pulse, -1 previous pulse
void ApplyWeightsMove(int direction) {
  int iPulse = aw_iPulseCurrent + direction;
  if (iPulse < 0 || iPulse >= int(aw_pulses.size())) {
    std::cout << " no " << (direction > 0 ? "next" : "previous") << " pulse with ADC[5] - pedestal > " << aw_threshold << std::endl;
    return;
  }
  aw_iPulseCurrent = iPulse;
  ApplyWeightsDraw(iPulse);
}


void ApplyWeights(std::string inputFile = "samples.root", std::string weightsFile = "weights_pulses_all.txt", float threshold = 50, int ieta_sel = 0, int iphi_sel = 0, int iEvent_sel = -1) {

  const int nSamples = aw_nSamples;
  aw_threshold = threshold;

  //---- read the weights, all ieta
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
    std::vector<float> weights;
    float value;
    while (iss >> value) weights.push_back(value);
    if (weights.size() == nSamples) aw_weights[ieta] = weights;
  }
  file.close();
  std::cout << " read weights for " << aw_weights.size() << " ieta from " << weightsFile << std::endl;
  if (ieta_sel != 0 && aw_weights.find(ieta_sel) == aw_weights.end()) {
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

  //---- find all the pulses of the selected crystals, one loop over the tree
  aw_pulses.clear();
  std::cout << " looking for the pulses ..." << std::endl;
  for (int iEvent = 0; iEvent < aw_tree->GetEntries(); iEvent++) {
    aw_tree->GetEntry(iEvent);
    for (int iDigi = 0; iDigi < aw_nDigis; iDigi++) {
      if (ieta_sel != 0 && aw_ieta[iDigi] != ieta_sel) continue;
      if (iphi_sel != 0 && aw_iphi[iDigi] != iphi_sel) continue;
      if (aw_weights.find(aw_ieta[iDigi]) == aw_weights.end()) continue;
      float pedestal = 0.5 * (aw_adc[iDigi*nSamples] + aw_adc[iDigi*nSamples + 1]);
      if (aw_adc[iDigi*nSamples + aw_iSampleMax] - pedestal <= threshold) continue;
      aw_pulses.push_back({aw_ieta[iDigi], aw_iphi[iDigi], iEvent, iDigi});
    }
  }

  //---- order: ieta (from -85), iphi (from 1), event
  std::sort(aw_pulses.begin(), aw_pulses.end(), [](const aw_pulse& a, const aw_pulse& b) {
    if (a.ieta != b.ieta) return a.ieta < b.ieta;
    if (a.iphi != b.iphi) return a.iphi < b.iphi;
    return a.iEvent < b.iEvent;
  });
  std::cout << " found " << aw_pulses.size() << " pulses with ADC[5] - pedestal > " << threshold << std::endl;
  if (aw_pulses.size() == 0) return;

  //---- first pulse: the first one, or the first one in event iEvent_sel
  aw_iPulseCurrent = 0;
  if (iEvent_sel >= 0) {
    aw_iPulseCurrent = -1;
    for (int iPulse = 0; iPulse < int(aw_pulses.size()); iPulse++) {
      if (aw_pulses[iPulse].iEvent == iEvent_sel) {
        aw_iPulseCurrent = iPulse;
        break;
      }
    }
    if (aw_iPulseCurrent < 0) {
      std::cout << " no pulse in event " << iEvent_sel << ", start from the first one" << std::endl;
      aw_iPulseCurrent = 0;
    }
  }

  aw_cc = new TCanvas ("cc", "", 800, 600);
  ApplyWeightsDraw(aw_iPulseCurrent);

  //---- control bar with the buttons, separate window as in the TMVA GUI
  TControlBar* bar = new TControlBar ("vertical", "ApplyWeights", 0, 0);
  bar->SetButtonWidth(200);
  bar->AddButton("+", "ApplyWeightsMove(1)",  "next pulse");
  bar->AddButton("-", "ApplyWeightsMove(-1)", "previous pulse");
  bar->Show();

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
