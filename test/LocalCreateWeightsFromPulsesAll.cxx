//
// Create weights for ECAL phase-2 barrel local reconstruction, one set per crystal,
// from the mean pulse of each crystal in the tree of SamplesTreeProducer
//   (clone of LocalCreateWeightsFromPulses.cxx, per crystal instead of per ieta ring)
//
//   pedestal = (ADC[0] + ADC[1]) / 2, subtracted to all the samples before anything else (also before the cut)
//   pulses: ADC[5] - pedestal > threshold, normalized as (ADC - pedestal) / (ADC[5] - pedestal)
//   mean pulse per crystal, only crystals with at least minPulses pulses (their number is printed)
//
//   W = (C^T N^-1 C)^-1 C^T N^-1
//
// with N the noise correlation matrix (unit matrix if useNoiseCorrelation = false)
// and C the (nSamples x nParams) design matrix:
//   column 0                 --> pulse shape (amplitude)
//   column 1 (doFitTime)     --> pulse shape derivative (amplitude * dt [ns]), from a cubic spline through the 16 samples
//   last col (doFitBaseline) --> 1 (baseline)
//
// Output file: amplitude weights, one line per crystal
//   # ieta iphi weights
//   ieta iphi w[0] ... w[15]
//
// Run:
//   root -l -q 'LocalCreateWeightsFromPulsesAll.cxx("samples.root", 50, "weights_pulses_crystals.txt")'
//   root -l -q 'LocalCreateWeightsFromPulsesAll.cxx("samples.root", 50, "weights_pulses_crystals.txt", 3, true, true, true)'
//

const int nSamples = 16;
const float samplingPeriod = 6.25;  // ns


//---- weights for one pulse shape, W = (C^T N^-1 C)^-1 C^T N^-1  (nParams x nSamples)
TMatrixD computeWeights(std::vector<float>& pulseShape, std::vector<float>& pulseShapeDerivative, bool doFitBaseline, bool doFitTime, bool useNoiseCorrelation, double& det) {

  //---- noise correlation G10 from 2023 TB, as in SimCalorimetry/EcalSimProducers/python/ecalDigiParameters_Ph2_cff.py
  float noiseCorrelation[nSamples] = {1.0000, 0.6986, 0.3782, 0.2452, 0.2900, 0.3532, 0.3911, 0.3896,
                                      0.3892, 0.3842, 0.3796, 0.3711, 0.3763, 0.3791, 0.3867, 0.3820};

  //---- derivative [1/ns] from a cubic spline through the samples
  double times[nSamples];
  double shape[nSamples];
  for (int iSample = 0; iSample < nSamples; iSample++) {
    times[iSample] = iSample * samplingPeriod;
    shape[iSample] = pulseShape[iSample];
  }
  TSpline3 spline ("spline", times, shape, nSamples);
  pulseShapeDerivative.clear();
  for (int iSample = 0; iSample < nSamples; iSample++) pulseShapeDerivative.push_back(spline.Derivative(times[iSample]));

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
  noiseInverse.Invert(&det);

  TMatrixD coefT (TMatrixD::kTransposed, coef);
  TMatrixD normal = coefT * noiseInverse * coef;
  normal.Invert(&det);
  return normal * coefT * noiseInverse;
}



void LocalCreateWeightsFromPulsesAll(std::string inputFile = "samples.root", float threshold = 50, std::string outputFile = "weights_pulses_crystals.txt", int minPulses = 3, bool doFitBaseline = true, bool doFitTime = true, bool useNoiseCorrelation = false) {

  const int iSampleMax = 5;   // sample at the peak (phase-2 simulation), counting from 0
  const int nEta = 85;        // ieta = -85 ... 85, no 0
  const int nPhi = 360;       // iphi = 1 ... 360

  std::cout << " Create phase-2 weights per crystal from: " << inputFile << ", ADC[5] - pedestal > " << threshold << std::endl;

  TFile* _file0 = TFile::Open(inputFile.c_str());
  TTree* tree = (TTree*) _file0->Get("SamplesTreeProducer/tree");

  int nDigis;
  static int      ieta[61200];
  static int      iphi[61200];
  static UShort_t adc[61200*16];

  tree->SetBranchStatus("*", 0);
  tree->SetBranchStatus("nDigis", 1);
  tree->SetBranchStatus("nSamplesTot", 1);
  tree->SetBranchStatus("ieta", 1);
  tree->SetBranchStatus("iphi", 1);
  tree->SetBranchStatus("adc", 1);

  tree->SetBranchAddress("nDigis", &nDigis);
  tree->SetBranchAddress("ieta",   ieta);
  tree->SetBranchAddress("iphi",   iphi);
  tree->SetBranchAddress("adc",    adc);

  //---- sum of the normalized pulses and number of pulses per crystal, index = (ieta + nEta) * (nPhi + 1) + iphi
  const int nCrystals = (2*nEta + 1) * (nPhi + 1);
  std::vector< std::vector<double> > sum_pulse (nCrystals, std::vector<double> (nSamples, 0.));
  std::vector<int> nPulses (nCrystals, 0);

  //---- loop over all the events
  for (int iEvent = 0; iEvent < tree->GetEntries(); iEvent++) {
    tree->GetEntry(iEvent);
    if (iEvent % 100 == 0) std::cout << " event " << iEvent << " / " << tree->GetEntries() << std::endl;

    for (int iDigi = 0; iDigi < nDigis; iDigi++) {
      if (abs(ieta[iDigi]) > nEta || iphi[iDigi] < 1 || iphi[iDigi] > nPhi) continue;

      //---- pedestal from the first two samples, subtracted before anything else
      float pedestal = 0.5 * (adc[iDigi*nSamples] + adc[iDigi*nSamples + 1]);
      //---- normalize to the sample at the peak
      float maximum = adc[iDigi*nSamples + iSampleMax] - pedestal;
      if (maximum <= 0 || maximum <= threshold) continue;

      int index = (ieta[iDigi] + nEta) * (nPhi + 1) + iphi[iDigi];
      nPulses[index]++;
      for (int iSample = 0; iSample < nSamples; iSample++) sum_pulse[index][iSample] += (adc[iDigi*nSamples + iSample] - pedestal) / maximum;
    }
  }

  //---- count the crystals with at least minPulses pulses
  int nCrystalsGood = 0;
  int nCrystalsAny = 0;
  for (int index = 0; index < nCrystals; index++) {
    if (nPulses[index] > 0) nCrystalsAny++;
    if (nPulses[index] >= minPulses) nCrystalsGood++;
  }
  std::cout << std::endl;
  std::cout << " crystals with at least 1 pulse with ADC[5] - pedestal > " << threshold << ": " << nCrystalsAny << std::endl;
  std::cout << " crystals with at least " << minPulses << " pulses with ADC[5] - pedestal > " << threshold << ": " << nCrystalsGood << " / " << 2*nEta*nPhi << std::endl;

  //---- weights per crystal
  std::ofstream out (outputFile.c_str());
  out << "# ieta iphi weights" << std::endl;
  int nWritten = 0;
  for (int iEta = -nEta; iEta <= nEta; iEta++) {
    if (iEta == 0) continue;
    for (int iPhi = 1; iPhi <= nPhi; iPhi++) {
      int index = (iEta + nEta) * (nPhi + 1) + iPhi;
      if (nPulses[index] < minPulses) continue;

      std::vector <float> pulseShape;
      for (int iSample = 0; iSample < nSamples; iSample++) pulseShape.push_back(sum_pulse[index][iSample] / nPulses[index]);
      std::vector <float> pulseShapeDerivative;
      double det = 0;
      TMatrixD weights = computeWeights(pulseShape, pulseShapeDerivative, doFitBaseline, doFitTime, useNoiseCorrelation, det);
      if (det == 0) {
        std::cout << " Error: singular matrix for ieta = " << iEta << ", iphi = " << iPhi << ", skipped" << std::endl;
        continue;
      }

      out << iEta << " " << iPhi;
      for (int iSample = 0; iSample < nSamples; iSample++) out << " " << std::setprecision(7) << weights(0, iSample);
      out << std::endl;
      nWritten++;
    }
  }
  out.close();
  std::cout << " written: " << outputFile << " (" << nWritten << " crystals)" << std::endl;

}
