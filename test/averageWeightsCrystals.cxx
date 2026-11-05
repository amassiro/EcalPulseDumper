//
// Average of the amplitude weights per crystal (from LocalCreateWeightsFromPulsesAll.cxx)
//   (clone of averageWeights.cxx, per crystal instead of per ieta ring)
//
//   input:  "# ieta iphi weights", then one line per crystal: ieta iphi w[0] ... w[15]
//   output: average_<input file name>, one line with the 16 averaged weights
//           canvas with the distributions of w[0] ... w[15] over the crystals,
//           saved in outputFolder as weights_distributions_crystals.png and .root
//
// Run:
//   root -l 'averageWeightsCrystals.cxx("weights_pulses_crystals.txt")'
//

void averageWeightsCrystals(std::string weightsFile = "weights_pulses_crystals.txt", std::string outputFolder = "pulses_output") {

  const int nSamples = 16;

  //---- read the weights, all crystals
  std::vector<double> average (nSamples, 0.);
  std::vector< std::vector<float> > all_weights;
  int nCrystals = 0;
  std::ifstream file (weightsFile.c_str());
  if (!file.is_open()) {
    std::cout << " Error: cannot open " << weightsFile << std::endl;
    return;
  }
  std::string line;
  while (std::getline(file, line)) {
    if (line.size() == 0 || line[0] == '#') continue;
    std::istringstream iss (line);
    int ieta, iphi;
    iss >> ieta >> iphi;
    std::vector<float> weights;
    float value;
    while (iss >> value) weights.push_back(value);
    if (weights.size() != nSamples) continue;
    for (int iSample = 0; iSample < nSamples; iSample++) average[iSample] += weights[iSample];
    all_weights.push_back(weights);
    nCrystals++;
  }
  file.close();
  if (nCrystals == 0) {
    std::cout << " Error: no weights found in " << weightsFile << std::endl;
    return;
  }
  for (int iSample = 0; iSample < nSamples; iSample++) average[iSample] /= nCrystals;

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

  std::cout << " average of " << nCrystals << " crystals:" << std::endl << "   ";
  for (int iSample = 0; iSample < nSamples; iSample++) std::cout << std::setprecision(5) << std::fixed << average[iSample] << " ";
  std::cout << std::endl;
  std::cout << " written: " << outputFile << std::endl;

  //---- distributions of the weights over the crystals, one histogram per sample
  gSystem->mkdir(outputFolder.c_str(), true);
  //---- bigger title, axis labels: the pads are small
  gStyle->SetTitleFontSize(0.09);
  TCanvas* cc = new TCanvas ("weights_distributions_crystals", "", 1600, 1200);
  cc->Divide(4,4);
  for (int iSample = 0; iSample < nSamples; iSample++) {
    float minimum = all_weights[0][iSample];
    float maximum = all_weights[0][iSample];
    for (auto& weights : all_weights) {
      if (weights[iSample] < minimum) minimum = weights[iSample];
      if (weights[iSample] > maximum) maximum = weights[iSample];
    }
    float margin = 0.1 * (maximum - minimum) + 1e-4;
    TH1F* h_weight = new TH1F (Form("h_weight_%d", iSample), Form("w[%d];w[%d];crystals", iSample, iSample), 30, minimum - margin, maximum + margin);
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
    h_weight->SetStats(0);
    h_weight->Draw();

    //---- own stat box: std dev with 2 significant digits, mean rounded to the same decimal place
    double stddev = h_weight->GetStdDev();
    int decimals = (stddev > 0) ? std::max(0, 1 - int(std::floor(std::log10(stddev)))) : 4;
    TPaveText* stat = new TPaveText (0.55, 0.70, 0.98, 0.92, "NDC");
    stat->SetFillColor(kWhite);
    stat->SetBorderSize(1);
    stat->SetTextAlign(12);
    stat->SetTextSize(0.07);
    stat->SetTextFont(42);
    stat->AddText(Form("Entries  %d", int(h_weight->GetEntries())));
    stat->AddText(Form("Mean  %.*f", decimals, h_weight->GetMean()));
    stat->AddText(Form("Std Dev  %.*f", decimals, stddev));
    stat->Draw();
  }
  cc->SaveAs((outputFolder + "/weights_distributions_crystals.png").c_str());
  cc->SaveAs((outputFolder + "/weights_distributions_crystals.root").c_str());

}
