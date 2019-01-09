#include "AveragePowerSpectrum.h"

ClassImp(AveragePowerSpectrum);





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Constructor
 *
 * For ROOT, Don't use this
 */
AveragePowerSpectrum::AveragePowerSpectrum(){
  // Don't use this.
  for(Int_t freqInd=0; freqInd<NUM_FREQS; freqInd++){  
    hRayleighs[freqInd] = NULL;
    // hRayleighFits[freqInd] = NULL;
  }
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Constructor
 * 
 * @param name is the name you want to give this average power spectrum (e.g. channel).
 * @param title is the title you want to give this average power spectrum, since we inherit from TNamed.
 */
AveragePowerSpectrum::AveragePowerSpectrum(TString name, TString title){


  deltaFMHz = 1e3/(NOMINAL_SAMPLING_DELTAT*NUM_SAMPLES);

  count=0;
  SetName(name);
  SetTitle(title);  

  // Default power spec histogram values.
  // maxNumOutliers = 5;
  // maxNumOutliers = 0;  

  for(Int_t freqInd=0; freqInd<NUM_FREQS; freqInd++){
    // psdOutliers[freqInd] = std::vector<Double_t>(0, 0);
    summedPowSpec[freqInd] = 0;

    TString histName = name + TString::Format("_%d", freqInd);
    TString histTitle = title + TString::Format(" Rayleigh Distribution for %4.2lf MHz bin",
						deltaFMHz*freqInd);
    histTitle += "; Amplitude (ADC/MHz); Events/bin";
    TH1D* hTemp = new TH1D(histName, histTitle,
			   NUM_AMPLITUDE_BINS,
			   INITIAL_MIN_AMPLITUDE,
			   INITIAL_MAX_AMPLITUDE);
    hTemp->SetDirectory(0);
    hTemp->Sumw2();

//     // This prepocessor variable is defined in the Makefile and is (I hope)
//     // querying the ROOT major version number. At the moment it asks whether the
//     // ROOT major version number is >= 6, hence the name.
//     // It works for me, for now.
// #ifdef IS_ROOT_6
    hTemp->SetCanExtend(TH1::kAllAxes);
// #else
//     hTemp->SetBit(TH1::kCanRebin);
// #endif

    hRayleighs[freqInd] = hTemp;
    // hDiffPhases[freqInd] = new TH1D(Form("hDiffPhases_%d", freqInd), "Differential Phases", NUM_AMPLITUDE_BINS, -TMath::Pi(), +TMath::Pi() );

    // hRayleighFits[freqInd] = NULL;


    summedDifferentialPhases[freqInd] = 0;
    summedDifferentialPhasesSq[freqInd] = 0;

    rayleighFitChiSquares[freqInd] = -1;
    rayleighFitChiSquaresRisingEdge[freqInd] = -1;
    rayleighFitChiSquaresRisingEdgeAndHalfFalling[freqInd] = -1;

    rayleighAmplitudes[freqInd] = -1;
    rayleighAmplitudesRisingEdge[freqInd] = -1;
    rayleighAmplitudesRisingEdgeAndHalfFalling[freqInd] = -1;

    rayleighNdf[freqInd] = -1;
    rayleighNdfRisingEdge[freqInd] = -1;
    rayleighNdfRisingEdgeAndHalfFalling[freqInd] = -1;

    for(int i=0; i<MAX_NUM_OUTLIERS; i++){
      outliers[freqInd][i] = -1;
    }
    numOutliers[freqInd] = 0;
  }
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Destructor
 */
AveragePowerSpectrum::~AveragePowerSpectrum(){
  deleteRayleighDistributions();
}




//---------------------------------------------------------------------------------------------------------
/**
 * @brief Generates the text that defines the Rayleigh function for TF1.
 */
TString AveragePowerSpectrum::getRayleighFunctionText(){
  return "([0]*x/([1]*[1]))*exp(-x*x/(2*[1]*[1]))";
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Creates a TF1 with to fit a Rayeligh distribution.
 *
 * @param name is the name for the TF1, use to identify the frequency bin.
 * @param xMin is the low edge of the fit range.
 * @param xMax is the high edge of the fit range.
 */
TF1* AveragePowerSpectrum::makeRayleighFunction(TString name, Double_t xMin, Double_t xMax){
  TF1* fRay = new TF1(name, getRayleighFunctionText(), xMin, xMax);
  return fRay;
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Converts the stored power spectrum for the events and stores it in internal memory
 *
 * @param gr is a TGraph containing the voltage/time waveform
 */
void AveragePowerSpectrum::getEventRayleighAmplitudes(TGraph* gr){


  // Original ANITA code used FancyFFTs
  // Double_t* ps = FancyFFTs::getPowerSpectrum(NUM_SAMPLES, gr->GetY(),
  // 					     NOMINAL_SAMPLING_DELTAT, FancyFFTs::kSum);

  Double_t thisPhase[NUM_FREQS];
  //  TGraph *gtempPow = FFTtools::makeRawPowerSpectrum(gr);
  TGraph *gtempPow = makeRawPowerSpectrum(gr, thisPhase);
  Int_t length = gtempPow->GetN();
  Double_t* ps = gtempPow->GetY();

  eventDiffPhases[0] = 0.;
  for(Int_t freqInd=0; freqInd < NUM_FREQS; freqInd++){
    // Double_t sqrtPSD = TMath::Sqrt(ps[freqInd]);

    ps[freqInd] /= length;
    Double_t sqrtPSD = TMath::Sqrt(ps[freqInd])/(deltaFMHz);

    eventRayleighAmplitudes[freqInd] = sqrtPSD;
    eventPowSpec[freqInd] = ps[freqInd];    
    if (freqInd>0) eventDiffPhases[freqInd] = FFTtools::wrap(thisPhase[freqInd]-thisPhase[freqInd-1], TMath::Pi()*2., 0.);
  }

  delete gtempPow;
  //  delete [] ps;
}



TGraph *AveragePowerSpectrum::makeRawPowerSpectrum(const TGraph *grWave, double phases[NUM_FREQS]) {

  double *oldY = grWave->GetY();
  double *oldX = grWave->GetX();
  double deltaT=oldX[1]-oldX[0];
  int length=grWave->GetN();
  FFTWComplex *theFFT=FFTtools::doFFT(length,oldY);

  int newLength=(length/2)+1;
  double *newY = new double [newLength];
  double *newX = new double [newLength];

  double deltaF=1/(deltaT*length);
  //    double fMax = 1/(2*deltaT);  // In GHz

  double tempF=0;
  double power=0;
  for(int i=0;i<newLength;i++) {
    //    float power=pow(getAbs(theFFT[i]),2);
    power       = theFFT[i].getAbsSq();
    phases[i]   = theFFT[i].getPhase();
    if(i>0 && i<newLength-1) power*=2; //account for symmetry
    newX[i]=tempF;
    newY[i]=power;
    tempF+=deltaF;
  }


  TGraph *grPower = new TGraph(newLength,newX,newY);
  delete [] theFFT;
  delete [] newY;
  delete [] newX;
  return grPower;

}



//---------------------------------------------------------------------------------------------------------
/**
 * @brief Function for the user to add the voltage/time waveform to all the stored averages.
 *
 * @param gr is a TGraph containing the voltage/time waveform.
 * @return the number of waveforms added to the AveragePowerSpectrum instance.
 *
 * This function stores the average power spectrum and rayleigh amplitudes in interal memory.
 */
size_t AveragePowerSpectrum::add(TGraph* gr){

  // Double_t* ps = FancyFFTs::getPowerSpectrum(numSamples, gr->GetY(), 1e-3*deltaT, FancyFFTs::kPowSpecDensity);
  // Double_t* ps = FancyFFTs::getPowerSpectrum(numSamples, gr->GetY(), deltaT, FancyFFTs::kPowSpecDensity);

  getEventRayleighAmplitudes(gr);
  
  for(Int_t freqInd=0; freqInd < NUM_FREQS; freqInd++){
    // Double_t sqrtPSD = TMath::Sqrt(ps[freqInd]);
    Double_t sqrtPSD = eventRayleighAmplitudes[freqInd];

    //    std::cout << sqrtPSD << std::endl;

    TH1D* h = hRayleighs[freqInd];
    Double_t histMaxVal = h->GetXaxis()->GetBinLowEdge(h->GetNbinsX()+1);
    if(sqrtPSD < histMaxVal){
      h->Fill(sqrtPSD);
    }
    else{
      outliers[freqInd][numOutliers[freqInd]] = sqrtPSD;
      numOutliers[freqInd]++;


      if(numOutliers[freqInd]>=MAX_NUM_OUTLIERS){

	while(numOutliers[freqInd] >= MAX_NUM_OUTLIERS){
	  histMaxVal*=2;
	    
	  Double_t outliersTemp[MAX_NUM_OUTLIERS] = {0};
	  Int_t numOutliersTemp = 0;
	  
	  for(int i=0; i<MAX_NUM_OUTLIERS; i++){
	    if(outliers[freqInd][i] < histMaxVal){
	      h->Fill(outliers[freqInd][i]);
	    }
	    else{
	      outliersTemp[numOutliersTemp] = outliers[freqInd][i];
	      numOutliersTemp++;
	    }
	  }

	  numOutliers[freqInd] = 0;
	  for(int i=0; i < numOutliersTemp; i++){
	    outliers[freqInd][i] = outliersTemp[i];
	    numOutliers[freqInd]++;
	  }
	  for(int i=numOutliersTemp; i < MAX_NUM_OUTLIERS; i++){
	    outliers[freqInd][i] = -1;
	  }
	}
      }
    }
  }

  for(Int_t freqInd=0; freqInd < NUM_FREQS; freqInd++){
    summedPowSpec[freqInd] += eventPowSpec[freqInd];
    summedDifferentialPhases[freqInd] += eventDiffPhases[freqInd];
    summedDifferentialPhasesSq[freqInd] += eventDiffPhases[freqInd]*eventDiffPhases[freqInd];

    // hDiffPhases[freqInd]->Fill(eventDiffPhases[freqInd]);

  }
  count++;
  
  return count;
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Rebins all the stored Rayleigh distributions.
 *
 * @param rebinFactor is the factor to reduce the number of bins in the histogram by.
 */
void AveragePowerSpectrum::rebinAllRayleighHistograms(Int_t rebinFactor){
  for(Int_t freqInd=0; freqInd < NUM_FREQS; freqInd++){
    TH1D* h = hRayleighs[freqInd];
    h->Rebin(rebinFactor);
  }
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Get pointer to rayleigh distribution of frequency closest to a particular frequency.
 *
 * @param freqMHz is the chosen frequency in MHz.
 * @return pointer to the stored Rayleigh distribution histogram.
 */
TH1D* AveragePowerSpectrum::getRayleighHistogramFromFrequencyMHz(Double_t freqMHz){

  Int_t bestFreqInd=0;
  Double_t bestFreqDiff = DBL_MAX;
  for(Int_t freqInd=0; freqInd < NUM_FREQS-1; freqInd++){
    Double_t freqDiff = TMath::Abs(freqInd*deltaFMHz - freqMHz);
    if(freqDiff < bestFreqDiff){
      bestFreqInd = freqInd;
      bestFreqDiff = freqDiff;
    }
  }
  return hRayleighs[bestFreqInd];
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Get Rayleigh distribution from frequency index.
 *
 * @param freqInd is the index of the frequency bin.
 * @return pointer to the stored Rayleigh distribution histogram.
 */
TH1D* AveragePowerSpectrum::getRayleighHistogram(Int_t freqInd){
  return hRayleighs[freqInd];
}

// TH1D *AveragePowerSpectrum::getDiffPhaseHistogram(Int_t freqInd){
//   return hDiffPhases[freqInd];
// }




//---------------------------------------------------------------------------------------------------------
/**
 * @brief Get Rayleigh fit to Rayleigh distribution from frequency index.
 *
 * @param freqInd is the index of the frequency bin.
 * @return pointer to the stored Rayleigh distribution histogram fit.
 */
TF1* AveragePowerSpectrum::getRayleighHistogramFit(Int_t freqInd){
  TH1D* h = hRayleighs[freqInd];
  TString funcName = TString::Format("fit_%s_%d",  GetName(), freqInd);
  TF1* f = (TF1*) h->FindObject(funcName);
  return f;
}






//---------------------------------------------------------------------------------------------------------
/**
 * @brief Combine all the 1D rayleigh distribution histograms to a 2D histogram.
 *
 * @return pointer to the new TH2D.
 */
TH2D* AveragePowerSpectrum::makeRayleigh2DHistogram(){

  Double_t maxVal = 0;
  for(Int_t freqInd=0; freqInd < NUM_FREQS; freqInd++){
    TH1D* h = getRayleighHistogram(freqInd);
    Double_t xMax = h->GetXaxis()->GetBinLowEdge(h->GetNbinsX()+1);
    if(xMax > maxVal){
      maxVal = xMax;
    }
  }


  TString name = TString::Format("h2D_%s", GetName());
  TString title = TString::Format("%s Rayleigh Distribution Summary", GetTitle());
  
  TH2D* h2 = new TH2D(name, title, NUM_FREQS, 0, NUM_FREQS*deltaFMHz,
		      NUM_AMPLITUDE_BINS, 0, maxVal);

  h2->GetXaxis()->SetTitle("Frequency (MHz)");
  h2->GetXaxis()->SetNoExponent(1);
  h2->GetYaxis()->SetTitle("Amplitude (ADC/MHz)");  
  h2->GetYaxis()->SetNoExponent(1);
  h2->Sumw2();
  
  for(Int_t freqInd=0; freqInd < NUM_FREQS; freqInd++){
    TH1D* h = getRayleighHistogram(freqInd);
    for(Int_t binx=1; binx<=h->GetNbinsX(); binx++){
      Double_t freqMHz = freqInd*deltaFMHz;
      Double_t amplitude = h->GetBinCenter(binx);
      Double_t weight = h->GetBinContent(binx);      
      
      h2->Fill(freqMHz, amplitude, weight);
    }
  }
  
  return h2;
  
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Fits all the Rayeligh distributions
 */
void AveragePowerSpectrum::fitAllRayleighHistograms(){
  for(Int_t freqInd=0; freqInd < NUM_FREQS; freqInd++){
    fitRayleighHistogram(freqInd);
  }
}






//---------------------------------------------------------------------------------------------------------
/**
 * @brief Fits all the Rayeligh distributions, but only up to the peak bin.
 */
void AveragePowerSpectrum::fitAllRayleighHistogramsRisingEdge(){
  for(Int_t freqInd=0; freqInd < NUM_FREQS; freqInd++){
    fitRayleighHistogramRisingEdge(freqInd);
  }
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Fits all the Rayeligh distributions, as far as the bin which is half the value of the peak bin.
 */
void AveragePowerSpectrum::fitAllRayleighHistogramsRisingEdgeAndHalfFallingEdge(){
  for(Int_t freqInd=0; freqInd < NUM_FREQS; freqInd++){
    fitRayleighHistogramRisingEdgeAndHalfFallingEdge(freqInd);
  }
}




//---------------------------------------------------------------------------------------------------------
/**
 * @brief Create a Rayleigh fit with a particular Rayleigh amplitude.
 *
 * @param freqInd is the bin of the frequency bin.
 * @param amplitude is the rayleigh amplitude.
 * @return a pointer to the new newly created Rayeligh fit.
 *
 * Limits on the fit come from the histogram.
 */
TF1* AveragePowerSpectrum::constructFitFromAmplitude(Int_t freqInd, Double_t amplitude){

  TH1D* h = getRayleighHistogram(freqInd);
  TString fitName = TString::Format("fit_%d_%lf", freqInd, amplitude);
  TF1* fit = makeRayleighFunction(fitName, 0, h->GetBinLowEdge(h->GetNbinsX()+1));
  fit->SetParameter(1, amplitude);
  fit->SetParameter(0, h->GetBinLowEdge(2)*h->Integral());  
  return fit;  
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Create a Rayleigh fit with stored Rayleigh amplitude.
 *
 * @param freqInd is the bin of the frequency bin.
 * @return a pointer to the new newly created Rayeligh fit.
 *
 * Limits on the fit come from the histogram.
 * The amplitude comes from the stored internal numbers.
 */
TF1* AveragePowerSpectrum::constructFitFromRayleighAmplitude(Int_t freqInd){
  return constructFitFromAmplitude(freqInd, rayleighAmplitudes[freqInd]);
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Create a Rayleigh fit with stored Rayleigh amplitude, fit up to the peak bin.
 *
 * @param freqInd is the bin of the frequency bin.
 * @return a pointer to the new newly created Rayeligh fit.
 *
 * Limits on the fit come from the histogram.
 * The amplitude comes from the stored internal numbers.
 */
TF1* AveragePowerSpectrum::constructFitFromRayleighAmplitudeRisingEdge(Int_t freqInd){
  return constructFitFromAmplitude(freqInd, rayleighAmplitudesRisingEdge[freqInd]);
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Create a Rayleigh fit with stored Rayleigh amplitude, fit past the peak to half the maximum value.
 *
 * @param freqInd is the bin of the frequency bin.
 * @return a pointer to the new newly created Rayeligh fit.
 *
 * Limits on the fit come from the histogram.
 * The amplitude comes from the stored internal numbers.
 */
TF1* AveragePowerSpectrum::constructFitFromRayleighAmplitudeRisingEdgeAndHalfFalling(Int_t freqInd){
  return constructFitFromAmplitude(freqInd, rayleighAmplitudesRisingEdgeAndHalfFalling[freqInd]);
}


//---------------------------------------------------------------------------------------------------------
/**
 * @brief Fit a Rayleigh histogram.
 *
 * @param freqInd is the bin of the frequency bin.
 */
void AveragePowerSpectrum::fitRayleighHistogram(Int_t freqInd){
  // std::cout << __PRETTY_FUNCTION__ << std::endl;
  
  TH1D* h = getRayleighHistogram(freqInd);
  xHigh[freqInd] = h->GetBinLowEdge(h->GetNbinsX()+1);

  fitRayleighHistogramOverRange(freqInd, h->GetBinLowEdge(1), xHigh[freqInd],
				rayleighAmplitudes,
				rayleighFitChiSquares,
				rayleighNdf,
				rayleighFitChiSquaresFullRange,
				rayleighNdfFullRange);
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Fit a Rayleigh histogram, up to the peak bin.
 *
 * @param freqInd is the bin of the frequency bin.
 */
void AveragePowerSpectrum::fitRayleighHistogramRisingEdge(Int_t freqInd){
  // std::cout << __PRETTY_FUNCTION__ << std::endl;

  TH1D* h = getRayleighHistogram(freqInd);
  //  Int_t peakBin = RootTools::getPeakBinOfHistogram(h);
  Int_t peakBin = h->GetMaximumBin();

  xHighRisingEdge[freqInd] = h->GetBinCenter(peakBin);
  fitRayleighHistogramOverRange(freqInd, h->GetBinLowEdge(1), xHighRisingEdge[freqInd],
				rayleighAmplitudesRisingEdge,
				rayleighFitChiSquaresRisingEdge,
				rayleighNdfRisingEdge,
				rayleighFitChiSquaresRisingEdgeFullRange,
				rayleighNdfRisingEdgeFullRange);
}






//---------------------------------------------------------------------------------------------------------
/**
 * @brief Fit a Rayleigh histogram, past the peak up to half the maximum value.
 *
 * @param freqInd is the bin of the frequency bin.
 */
void AveragePowerSpectrum::fitRayleighHistogramRisingEdgeAndHalfFallingEdge(Int_t freqInd){
  // std::cout << __PRETTY_FUNCTION__ << std::endl;

  TH1D* h = getRayleighHistogram(freqInd);
  Double_t xHighTemp = -1;
  Int_t peakBin = h->GetMaximumBin(); //RootTools::getPeakBinOfHistogram(h);
  Double_t peakVal = h->GetBinContent(peakBin);

  for(Int_t binx=peakBin; binx<=h->GetNbinsX(); binx++){
    Double_t binVal = h->GetBinContent(binx);
    if(binVal < peakVal*0.5){
      xHighTemp = h->GetBinCenter(binx);
      // std::cout << fName << "\t" << freqInd << "\t" << binx << "\t" << xHigh << std::endl;
      break;
    }
  }

  xHighRisingEdgeAndHalfFalling[freqInd] = xHighTemp;
  fitRayleighHistogramOverRange(freqInd, h->GetBinLowEdge(1), xHighTemp,
				rayleighAmplitudesRisingEdgeAndHalfFalling,
				rayleighFitChiSquaresRisingEdgeAndHalfFalling,
				rayleighNdfRisingEdgeAndHalfFalling,
				rayleighFitChiSquaresRisingEdgeAndHalfFallingFullRange,
				rayleighNdfRisingEdgeAndHalfFallingFullRange);
  // fitRayleighHistogramOverRange(freqInd, h->GetBinLowEdge(1), h->GetBinLowEdge(3));  
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Worker function called by all the other fitting functions.
 *
 * @param freqInd is the bin of the frequency bin.
 * @param xLowVal is the low range of the fit.
 * @param xHighVal is the high range of the fit.
 * @param rAmplitudes is a pointer to stored amptlitude array.
 * @param rChiSquares is a pointer to stored chi square array.
 * @param rNdf is a pointer to stored NDF array.
 * @param rChiSquaresFullRange is a pointer to stored chi square, over the whole range even if this is larger than the fitted range.
 * @param rNdfFullRange is a pointer to stored NDF, over the whole range even if this is larger than the fitted range.
 */
void AveragePowerSpectrum::fitRayleighHistogramOverRange(Int_t freqInd, Double_t xLowVal, Double_t xHighVal,
							 Double_t* rAmplitudes,
							 Double_t* rChiSquares,
							 Int_t* rNdf,
							 Double_t* rChiSquaresFullRange,
							 Int_t* rNdfFullRange){


  TH1D* h = getRayleighHistogram(freqInd);

  if(h->Integral() > 0){ // Fit will fail with empty histogram

    TString fitName = TString::Format("fit_%s_%d", GetName(), freqInd );
    Double_t mean = h->GetMean();
    // Mean of rayleigh distribution = sigma* (pi/2)^{0.5}
    Double_t sigGuess = mean / TMath::Sqrt(0.5*TMath::Pi());
    Double_t histArea = h->Integral() * h->GetBinLowEdge(2);
    
    //    TF1* fit;
// #pragma omp critical (fit)
//     {
      
      TF1 *fit = makeRayleighFunction(fitName, xLowVal, xHighVal);
      
      fit->SetParameter(1, sigGuess);
      
      // Can't not set upper bound if setting lower bound..
      // Therefore set the upper bound to be insanely high.
      // fit->SetParLimits(1, 0, h->GetMean()*1e9); 
      
      fit->FixParameter(0, histArea);
      
      h->Fit(fit, "RQ0");
      //    }
    // h->GetFunction(fitName.Data())->ResetBit(TF1::kNotDraw);
    rAmplitudes[freqInd] = h->GetFunction(fitName.Data())->GetParameter(1);
    rChiSquares[freqInd] = h->GetFunction(fitName.Data())->GetChisquare();
    rNdf[freqInd] = h->GetFunction(fitName.Data())->GetNDF(); 
    
    rChiSquaresFullRange[freqInd] = 0;
    rNdfFullRange[freqInd] = 0;
    // for(int binx=1; binx<=h->GetNbinsX(); binx++){
    //   Double_t binVal = h->GetBinContent(binx);
    //   Double_t binError = h->GetBinError(binx);
    //   Double_t fitVal = fit->Eval(h->GetBinCenter(binx));
    //   if(binError > 0){
    // 	Double_t chi = (binVal - fitVal)/binError;
    // 	rChiSquaresFullRange[freqInd] += chi*chi;
    // 	rNdfFullRange[freqInd]++;
    //   }
    // }
    
    delete fit;
  }

  
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Deletes all the Rayleigh histograms
 */
void AveragePowerSpectrum::deleteRayleighDistributions(){
  for(int freqInd=0; freqInd < NUM_FREQS; freqInd++){
    if(hRayleighs[freqInd]!=NULL){
      delete hRayleighs[freqInd];
      hRayleighs[freqInd] = NULL;
    }
  }
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Creates a TGraph of the average power spectrum
 */
TGraph* AveragePowerSpectrum::makeAvePowSpecTGraph(){

  TString name = TString::Format("gr_%s", GetName());
  TString title = TString::Format("%s", GetTitle());
  
  std::vector<Double_t> avePowSpec(NUM_FREQS);
  
  for(Int_t freqInd=0; freqInd < NUM_FREQS; freqInd++){
    if(count > 0){
      avePowSpec[freqInd] = summedPowSpec[freqInd]/count;
    }
  }
  
  // Double_t termResOhms = 50;
  for(Int_t freqInd=0; freqInd < NUM_FREQS; freqInd++){
    avePowSpec[freqInd]/=(deltaFMHz);
  }
  
  Double_t freqArray[NUM_FREQS];
  for(int freqInd=0; freqInd < NUM_FREQS; freqInd++){
    freqArray[freqInd] = deltaFMHz*freqInd;
  }
  
  TGraph* gr = new TGraph(NUM_FREQS, freqArray, &avePowSpec[0]);
  gr->SetName(name);
  gr->SetTitle(title);
  return gr;
}






//---------------------------------------------------------------------------------------------------------
/**
 * @brief Creates a TGraph of the average power spectrum with a dB scale.
 */
TGraph* AveragePowerSpectrum::makeAvePowSpecTGraph_dB(){
  
  TGraph* gr = makeAvePowSpecTGraph();
  TString name = TString::Format("%s_dB", gr->GetName());

  for(Int_t freqInd=0; freqInd < NUM_FREQS; freqInd++){
    Double_t y = gr->GetY()[freqInd];
    // gr->GetY()[freqInd] = 10*TMath::Log10(y);
    gr->GetY()[freqInd] = 10*TMath::Log10(y);
  }
  return gr;
}
