//
// Average of the amplitude weights per ieta ring (from LocalCreateWeightsFromPulses.cxx)
//
//   input:  "# ieta weights", then one line per ieta: ieta w[0] ... w[15]
//   output: average_<input file name>, one line with the 16 averaged weights
//           canvas with the distributions of w[0] ... w[15] over the ieta rings,
//           saved in outputFolder as weights_distributions.png and .root
//
// Run:
//   root -l 'averageWeights.cxx("weights_pulses_all.txt")'
//

void averageWeights(std::string weightsFile = "weights_pulses_all.txt", std::string outputFolder = "pulses_output") {

  const int nSamples = 16;

  //---- read the weights, all ieta
  std::vector<double> average (nSamples, 0.);
  std::vector< std::vector<float> > all_weights;
  int nRings = 0;
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
    if (weights.size() != nSamples) continue;
    for (int iSample = 0; iSample < nSamples; iSample++) average[iSample] += weights[iSample];
    all_weights.push_back(weights);
    nRings++;
  }
  file.close();
  if (nRings == 0) {
    std::cout << " Error: no weights found in " << weightsFile << std::endl;
    return;
  }
  for (int iSample = 0; iSample < nSamples; iSample++) average[iSample] /= nRings;

  //---- output file: average_ + input file name, in the same folder
  std::string folder = "";
  std::string name = weightsFile;
  size_t slash = weightsFile.find_last_of('/');
  if (slash != std::string::npos) {
    folder = weightsFile.substr(0, slash + 1);
    name   = weightsFile.substr(slash + 1);
  }
  std::string outputFile = folder + "average_" + name;

  std::ofstream out (outputFile.c_str());
  for (int iSample = 0; iSample < nSamples; iSample++) out << std::setprecision(7) << average[iSample] << (iSample < nSamples - 1 ? " " : "");
  out << std::endl;
  out.close();

  std::cout << " average of " << nRings << " ieta rings:" << std::endl << "   ";
  for (int iSample = 0; iSample < nSamples; iSample++) std::cout << std::setprecision(5) << std::fixed << average[iSample] << " ";
  std::cout << std::endl;
  std::cout << " written: " << outputFile << std::endl;

  //---- distributions of the weights over the ieta rings, one histogram per sample
  gSystem->mkdir(outputFolder.c_str(), true);
  //---- bigger stat box, title, axis labels: the pads are small
  gStyle->SetOptStat(1110);
  gStyle->SetStatW(0.40);
  gStyle->SetStatH(0.25);
  gStyle->SetStatFontSize(0.07);
  //---- mean and std dev with 4 significant digits, trailing zeros kept
  gStyle->SetStatFormat("#.4g");
  gStyle->SetTitleFontSize(0.09);
  TCanvas* cc = new TCanvas ("weights_distributions", "", 1600, 1200);
  cc->Divide(4,4);
  for (int iSample = 0; iSample < nSamples; iSample++) {
    float minimum = all_weights[0][iSample];
    float maximum = all_weights[0][iSample];
    for (auto& weights : all_weights) {
      if (weights[iSample] < minimum) minimum = weights[iSample];
      if (weights[iSample] > maximum) maximum = weights[iSample];
    }
    float margin = 0.1 * (maximum - minimum) + 1e-4;
    TH1F* h_weight = new TH1F (Form("h_weight_%d", iSample), Form("w[%d];w[%d];ieta rings", iSample, iSample), 30, minimum - margin, maximum + margin);
    for (auto& weights : all_weights) h_weight->Fill(weights[iSample]);
    h_weight->SetMaximum(1.6 * h_weight->GetMaximum());
    h_weight->SetLineColor(kRed);
    h_weight->SetLineWidth(2);
    h_weight->GetXaxis()->SetLabelSize(0.06);
    h_weight->GetXaxis()->SetTitleSize(0.07);
    h_weight->GetXaxis()->SetNdivisions(505);
    h_weight->GetYaxis()->SetLabelSize(0.06);
    h_weight->GetYaxis()->SetTitleSize(0.07);
    h_weight->GetYaxis()->SetTitleOffset(0.8);
    cc->cd(iSample + 1);
    gPad->SetLeftMargin(0.14);
    gPad->SetBottomMargin(0.15);
    h_weight->Draw();
  }
  cc->SaveAs((outputFolder + "/weights_distributions.png").c_str());
  cc->SaveAs((outputFolder + "/weights_distributions.root").c_str());

}
