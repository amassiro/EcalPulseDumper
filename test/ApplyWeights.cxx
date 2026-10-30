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
// Run:
//   root -l 'ApplyWeights.cxx("samples.root", "weights_TB2025_v1.txt")'
//

void ApplyWeights(std::string inputFile = "samples.root", std::string weightsFile = "weights_Phase2.txt", float minAmplitude = 50) {

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
