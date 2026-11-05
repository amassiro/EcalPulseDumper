//
// Create weights for ECAL phase-2 barrel local reconstruction from the mean pulses measured per ieta
//   - input: text file from drawPulseAllEta.cxx (pulses_output/pulses.txt)
//            one line per ieta: ieta, ADC[0], ... , ADC[15]
//            with ADC = (ADC - pedestal) / (ADC[5] - pedestal) mean, pedestal = (ADC[0] + ADC[1]) / 2
//   - 16 samples at 160 MHz (6.25 ns)
//
//   W = (C^T N^-1 C)^-1 C^T N^-1
//
// with N the noise correlation matrix (unit matrix if useNoiseCorrelation = false)
// and C the (nSamples x nParams) design matrix:
//   column 0                 --> pulse shape (amplitude)
//   column 1 (doFitTime)     --> pulse shape derivative (amplitude * dt [ns]), from a cubic spline through the 16 samples
//                                A f(t) + A dt f'(t) ~ A f(t + dt) --> dt > 0 means the pulse is earlier
//   last col (doFitBaseline) --> 1 (baseline)
//
// ieta_sel = 0  --> mean of the pulses of all the ieta in the file
// ieta_sel != 0 --> pulse of that ieta
//
// Output: weights in a text file, one line per parameter (amplitude, [amplitude*dt], [baseline]),
//         same format as EcalWeights/LocalCreateWeightsPhase2.cxx (input of ApplyWeights.cxx)
//
// Run:
//   root -l -q 'LocalCreateWeightsFromPulses.cxx("pulses_output/pulses.txt", 0, "weights_pulses_all.txt")'
//   root -l -q 'LocalCreateWeightsFromPulses.cxx("pulses_output/pulses.txt", -54, "weights_pulses_ieta_m54.txt", true, true, true)'
//


void LocalCreateWeightsFromPulses(std::string pulsesFile = "pulses.txt", int ieta_sel = 0, std::string outputFile = "weights_pulses.txt", bool doFitBaseline = true, bool doFitTime = true, bool useNoiseCorrelation = false) {

  std::cout << " Create phase-2 weights from the mean pulses: " << pulsesFile << ", ieta = " << ieta_sel << " (0 = all)" << std::endl;

  const int nSamples = 16;
  const float samplingPeriod = 6.25;  // ns

  //---- noise correlation G10 from 2023 TB, as in SimCalorimetry/EcalSimProducers/python/ecalDigiParameters_Ph2_cff.py
  float noiseCorrelation[nSamples] = {1.0000, 0.6986, 0.3782, 0.2452, 0.2900, 0.3532, 0.3911, 0.3896,
                                      0.3892, 0.3842, 0.3796, 0.3711, 0.3763, 0.3791, 0.3867, 0.3820};

  //---- read the pulses: ieta  ADC[0] ... ADC[15]
  std::vector <float> pulseShape (nSamples, 0.);
  int nPulses = 0;
  std::ifstream file (pulsesFile.c_str());
  if (!file.is_open()) {
    std::cout << " Error: cannot open " << pulsesFile << std::endl;
    return;
  }
  std::string line;
  while (std::getline(file, line)) {
    if (line.size() == 0 || line[0] == '#') continue;
    std::istringstream iss (line);
    int ieta;
    iss >> ieta;
    if (ieta_sel != 0 && ieta != ieta_sel) continue;
    std::vector <float> values;
    float value;
    while (iss >> value) values.push_back(value);
    if (int(values.size()) != nSamples) continue;
    //---- skip empty rings (no pulses above threshold)
    if (values[5] == 0) continue;
    for (int iSample = 0; iSample < nSamples; iSample++) pulseShape[iSample] += values[iSample];
    nPulses++;
  }
  file.close();
  if (nPulses == 0) {
    std::cout << " Error: no pulse found for ieta = " << ieta_sel << std::endl;
    return;
  }
  for (int iSample = 0; iSample < nSamples; iSample++) pulseShape[iSample] /= nPulses;
  std::cout << " used " << nPulses << " ieta" << std::endl;

  //---- derivative [1/ns] from a cubic spline through the samples
  double times[nSamples];
  double shape[nSamples];
  for (int iSample = 0; iSample < nSamples; iSample++) {
    times[iSample] = iSample * samplingPeriod;
    shape[iSample] = pulseShape[iSample];
  }
  TSpline3* spline = new TSpline3 ("spline", times, shape, nSamples);
  std::vector <float> pulseShapeDerivative;
  for (int iSample = 0; iSample < nSamples; iSample++) pulseShapeDerivative.push_back(spline->Derivative(times[iSample]));

  int nParams = 1 + int(doFitBaseline) + int(doFitTime);
  int iParamTime     = doFitTime     ? 1           : -1;
  int iParamBaseline = doFitBaseline ? nParams - 1 : -1;

  //---- design matrix C (nSamples x nParams)
  TMatrixD coef (nSamples, nParams);
  coef.Zero();
  for (int iSample = 0; iSample < nSamples; iSample++) {
    coef(iSample, 0) = pulseShape[iSample];
    if (doFitTime)     coef(iSample, iParamTime) = pulseShapeDerivative[iSample];
    if (doFitBaseline) coef(iSample, iParamBaseline) = 1.;
  }

  //---- inverse noise matrix
  TMatrixD noiseInverse (nSamples, nSamples);
  for (int iRow = 0; iRow < nSamples; iRow++) {
    for (int iColumn = 0; iColumn < nSamples; iColumn++) {
      if (useNoiseCorrelation) noiseInverse(iRow, iColumn) = noiseCorrelation[abs(iRow - iColumn)];
      else                     noiseInverse(iRow, iColumn) = (iRow == iColumn) ? 1. : 0.;
    }
  }
  double det = 0;
  noiseInverse.Invert(&det);

  //---- W = (C^T N^-1 C)^-1 C^T N^-1   (nParams x nSamples)
  TMatrixD coefT (TMatrixD::kTransposed, coef);
  TMatrixD normal = coefT * noiseInverse * coef;
  normal.Invert(&det);
  if (det == 0) {
    std::cout << " Error: singular matrix, check the pulse shape" << std::endl;
    return;
  }
  TMatrixD weights = normal * coefT * noiseInverse;

  //---- print
  std::cout << " noise correlation = " << useNoiseCorrelation << std::endl;
  std::cout << " parameters: 0 = amplitude";
  if (doFitTime)     std::cout << ", " << iParamTime << " = amplitude * dt [ns]";
  if (doFitBaseline) std::cout << ", " << iParamBaseline << " = baseline";
  std::cout << std::endl;

  std::cout << std::endl << " pulse shape: " << std::endl << "   ";
  for (int iSample = 0; iSample < nSamples; iSample++) std::cout << std::setw(8) << std::setprecision(4) << std::fixed << pulseShape[iSample] << " ";
  std::cout << std::endl;
  std::cout << " derivative [1/ns]: " << std::endl << "   ";
  for (int iSample = 0; iSample < nSamples; iSample++) std::cout << std::setw(8) << std::setprecision(4) << std::fixed << pulseShapeDerivative[iSample] << " ";
  std::cout << std::endl;

  std::cout << std::endl << " weights: " << std::endl;
  for (int iParam = 0; iParam < nParams; iParam++) {
    std::cout << "   [" << iParam << "] ";
    for (int iSample = 0; iSample < nSamples; iSample++) {
      std::cout << std::setw(9) << std::setprecision(5) << std::fixed << weights(iParam, iSample) << " ";
    }
    std::cout << std::endl;
  }

  //---- closure: weights applied to the pulse shape (+ constant baseline)
  std::cout << std::endl << " closure (weights x (pulse shape + 100)): " << std::endl;
  for (int iParam = 0; iParam < nParams; iParam++) {
    double sum = 0;
    for (int iSample = 0; iSample < nSamples; iSample++) sum += weights(iParam, iSample) * (pulseShape[iSample] + 100.);
    std::cout << "   [" << iParam << "] = " << sum << std::endl;
  }

  //---- cfg format, to be used in EcalUncalibRecHitPhase2WeightsProducer
  //---- NB: in the producer t0 = sum(timeWeights * adc), not divided by the amplitude
  std::cout << std::endl << " ampWeights = cms.vdouble(";
  for (int iSample = 0; iSample < nSamples; iSample++) std::cout << std::setprecision(6) << weights(0, iSample) << (iSample < nSamples - 1 ? ", " : ")");
  std::cout << std::endl;
  if (doFitTime) {
    std::cout << " timeWeights (amplitude * dt [ns]) = cms.vdouble(";
    for (int iSample = 0; iSample < nSamples; iSample++) std::cout << std::setprecision(6) << weights(iParamTime, iSample) << (iSample < nSamples - 1 ? ", " : ")");
    std::cout << std::endl;
  }

  //---- write weights to file
  std::ofstream out (outputFile.c_str());
  out << "# pulses " << pulsesFile << " ieta " << ieta_sel << " (0 = all) noiseCorrelation " << useNoiseCorrelation << std::endl;
  out << "# rows: amplitude";
  if (doFitTime)     out << ", amplitude*dt[ns]";
  if (doFitBaseline) out << ", baseline";
  out << std::endl;
  for (int iParam = 0; iParam < nParams; iParam++) {
    for (int iSample = 0; iSample < nSamples; iSample++) out << std::setprecision(7) << weights(iParam, iSample) << " ";
    out << std::endl;
  }
  out.close();
  std::cout << std::endl << " written: " << outputFile << std::endl;

}
