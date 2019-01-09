/* -*- C++ -*-.***************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk
 Description: 
 A Class to manage power spectra, probably needs to do things like rolling averages.
 
 Edited by Linda Cremonesi (l.cremonesi@ucl.ac.uk) to be used in protoDUNE
 *************************************************************************************************************** */

#ifndef AVERAGEPOWERSPECTRUM_H
#define AVERAGEPOWERSPECTRUM_H

#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TNamed.h"

#include "FFTtools.h"

/* #include "FancyFFTs.h" */
/* #include "RootTools.h" */
/* #include "CrossCorrelator.h" */

#define NUM_AMPLITUDE_BINS 64
#define INITIAL_MAX_AMPLITUDE 1
#define INITIAL_MIN_AMPLITUDE 0

#define MAX_NUM_OUTLIERS 10

#define NOMINAL_SAMPLING_DELTAT 500 // in ns

#define NUM_SAMPLES 6000
#define NUM_FREQS ((NUM_SAMPLES/2)+1)


/**
 * @class AveragePowerSpectrum
 * @brief Takes in waveforms and averages averages their power spectra
 */
class AveragePowerSpectrum : public TNamed {

 public:
  
  AveragePowerSpectrum();
  AveragePowerSpectrum(TString name, TString title);
  ~AveragePowerSpectrum();

  size_t add(TGraph* gr);
  
  TGraph* makeAvePowSpecTGraph();
  TGraph* makeAvePowSpecTGraph_dB();
  
  void deleteRayleighDistributions();
  void rebinAllRayleighHistograms(Int_t rebinFactor);

  
  void fitRayleighHistogramOverRange(Int_t freqInd, Double_t xLowVal, Double_t xHighVal,
				     Double_t* rAmplitudes,
				     Double_t* rChiSquares,
				     Int_t* rNdf,
				     Double_t* rChiSquaresFullRange,
				     Int_t* rNdfFullRange);


  void fitRayleighHistogram(Int_t freqInd);
  void fitRayleighHistogramRisingEdge(Int_t freqInd);
  void fitRayleighHistogramRisingEdgeAndHalfFallingEdge(Int_t freqInd);

  void fitAllRayleighHistograms();
  void fitAllRayleighHistogramsRisingEdge();
  void fitAllRayleighHistogramsRisingEdgeAndHalfFallingEdge();

  TH1D* getRayleighHistogram(Int_t freqInd);
  TH1D* getRayleighHistogramFromFrequencyMHz(Double_t freqMHz);
  TF1* getRayleighHistogramFit(Int_t freqInd);  

  TH2D* makeRayleigh2DHistogram();
  

  void getEventRayleighAmplitudes(TGraph* gr);
  
  static TF1* makeRayleighFunction(TString name, Double_t xMin, Double_t xMax);
  static TString getRayleighFunctionText();
  TF1* constructFitFromAmplitude(Int_t freqInd, Double_t amplitude);
  TF1* constructFitFromRayleighAmplitude(Int_t freqInd);
  TF1* constructFitFromRayleighAmplitudeRisingEdge(Int_t freqInd);
  TF1* constructFitFromRayleighAmplitudeRisingEdgeAndHalfFalling(Int_t freqInd);  

  TGraph *makeRawPowerSpectrum(const TGraph *grWave, double phases[NUM_FREQS]);

  Double_t deltaFMHz; //!< Difference between frequency bins in the power spectrums.
  Double_t summedPowSpec[NUM_FREQS]; //!< Sum of all power spectrum.
  TH1D* hRayleighs[NUM_FREQS]; //!< Histograms for Rayleigh distributions.

  // /////// TEMP
  // TH1D *hDiffPhases[NUM_FREQS];
  // TH1D *getDiffPhaseHistogram(Int_t freqInd);
  // //////// TEMP

  Int_t count; //!< Number of waveforms that have been added to the AveragePowerSpectrum.
  
  Double_t summedDifferentialPhases[NUM_FREQS]; //!< Sum of all differential phases
  Double_t summedDifferentialPhasesSq[NUM_FREQS]; //!< Sum of all differential phases

  Double_t rayleighFitChiSquares[NUM_FREQS]; //!< Chi squares of the Rayleigh fits.
  Double_t rayleighFitChiSquaresRisingEdge[NUM_FREQS]; //!< Chi squares of the fit to the rising egde of the Rayleigh distributions.
  Double_t rayleighFitChiSquaresRisingEdgeAndHalfFalling[NUM_FREQS]; //!< Chi squares of the fit to the leading and half the falling edge of the Rayleigh distributions.

  Double_t rayleighAmplitudes[NUM_FREQS]; //!< Fitted Rayleigh amplitudes
  Double_t rayleighAmplitudesRisingEdge[NUM_FREQS]; //!< Amplitudes from fits to the leading edge of the Rayleigh distributions.
  Double_t rayleighAmplitudesRisingEdgeAndHalfFalling[NUM_FREQS]; //!< Amplitudes from fits to the leading and half the falling edge of the Rayleigh distributions.

  Int_t rayleighNdf[NUM_FREQS]; //!< NDFs from Rayleigh fits
  Int_t rayleighNdfRisingEdge[NUM_FREQS]; //!< NDFs from fits to the leading edge of the Rayleigh distributions.
  Int_t rayleighNdfRisingEdgeAndHalfFalling[NUM_FREQS]; //!< NDFs from fits to the leading and half the falling edge of the Rayleigh distributions.

  Double_t rayleighFitChiSquaresFullRange[NUM_FREQS]; //!< Chi squares of the Rayleigh fits extended over the full range.
  Double_t rayleighFitChiSquaresRisingEdgeFullRange[NUM_FREQS]; //!< Chi squares of the fit to the rising egde of the Rayleigh distributions extended over the full range.
  Double_t rayleighFitChiSquaresRisingEdgeAndHalfFallingFullRange[NUM_FREQS]; //!< Amplitudes from fits to the leading and half the falling edge of the Rayleigh distributions extended over the full range.

  Int_t rayleighNdfFullRange[NUM_FREQS]; //!< NDFs from Rayleigh fits extended over the full range.
  Int_t rayleighNdfRisingEdgeFullRange[NUM_FREQS]; //!< NDFs from fits to the leading edge of the Rayleigh distributions extended over the full range.
  Int_t rayleighNdfRisingEdgeAndHalfFallingFullRange[NUM_FREQS]; //!< NDFs from fits to the leading and half the falling edge of the Rayleigh distributions extended over the full range.

  Double_t xHigh[NUM_FREQS]; //!< High end of fitted range of Rayleigh histograms.
  Double_t xHighRisingEdge[NUM_FREQS]; //!< High end of fitted range of Rayleigh histograms up to peak bin.
  Double_t xHighRisingEdgeAndHalfFalling[NUM_FREQS]; //!< High end of fitted range of Rayleigh histograms up to half the maximum value past the peak bin.

  Double_t outliers[NUM_FREQS][MAX_NUM_OUTLIERS]; //!< Used to store values that are much greater than the current Rayliegh histogram content.
  Int_t numOutliers[NUM_FREQS]; //!< Count of number of outliers in outliers array.

  Double_t eventRayleighAmplitudes[NUM_FREQS]; //!< Rayleigh amplitudes for the last added event.
  Double_t eventPowSpec[NUM_FREQS]; //!< Power spectra for the last added event.
  Double_t eventDiffPhases[NUM_FREQS]; //!< Differential phases for the last added event. 

  ClassDef(AveragePowerSpectrum, 11); //!< ROOT's magic I/O macro
};

  



#endif







