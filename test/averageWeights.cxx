//
// Average of the amplitude weights per ieta ring (from LocalCreateWeightsFromPulses.cxx)
//
//   input:  "# ieta weights", then one line per ieta: ieta w[0] ... w[15]
//   output: average_<input file name>, one line with the 16 averaged weights
//
// Run:
//   root -l -q 'averageWeights.cxx("weights_pulses_all.txt")'
//

void averageWeights(std::string weightsFile = "weights_pulses_all.txt") {

  const int nSamples = 16;

  //---- read the weights, all ieta
  std::vector<double> average (nSamples, 0.);
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

}
