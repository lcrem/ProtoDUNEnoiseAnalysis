#include "AveragePowerSpectrum.h"
#include "canvas/Utilities/InputTag.h"
#include "gallery/Event.h"
#include "lardataobj/RawData/RawDigit.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH2.h"
#include "TColor.h"
#include "TMath.h"
#include "TStyle.h"
#include "TSystem.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>

using namespace art;
using namespace std;
using namespace std::chrono;


int main(int argc, char *argv[]){

  Int_t run;
  Int_t subrun;
  Int_t channel;
  string inFileName;
  string outFileName;
  Int_t nevents=1e9;
  
  if((argc<2)){
    std::cerr << "Usage 1: " << argv[0] << " [irun] [subrun] [channel] (inlist) (outfile) (maxEvents)" << std::endl;
    return 1;
  } else {
    run = atoi(argv[1]);
    subrun = atoi(argv[2]);
    channel = atoi(argv[3]);
    if (argc>4) inFileName += argv[4];
    else  inFileName +=  "../filenameList.txt";
    if (argc>5) outFileName += argv[5];
    else outFileName = "tempRayleighTree.root";
    if (argc>6)  nevents = atoi(argv[6]);

  }

  cout << " Processing " << inFileName << endl;

  double x[NUM_SAMPLES];
  
  double dt = NOMINAL_SAMPLING_DELTAT;
  double deltaf = 2./NUM_SAMPLES;

  for (int i=0; i<NUM_SAMPLES; i++) x[i] = i*dt;

  TGraph *sampleGraph;
  TGraph *sampleFFT;
  bool doneGraph=false;

  InputTag const daq_tag{ "tpcrawdecoder", "daq", "DecoderandReco" }; 

  // Create a vector of length 1, containing the given filename.                                                                                   
  double samples[NUM_SAMPLES];

  ifstream FileList(inFileName.c_str());
  string filename;
  vector<string> filenames;
  while(FileList >> filename) {
    filenames.push_back(filename);
  }
  
  std::cout << "Have got a vector of filenames with length " << filenames.size() << "\n";

  int iev=0;
  double mean=0;
  int countsamples=0;


  
  Double_t powX[NUM_FREQS];
  Double_t rAmplitudes[NUM_FREQS];
  Double_t avgPowSpecY[NUM_FREQS];
  Double_t avgDiffPhases[NUM_FREQS];
  Double_t rmsDiffPhases[NUM_FREQS];
  Double_t rChiSquares[NUM_FREQS];
  Double_t rChiSquaresOverNdf[NUM_FREQS];
  Int_t rNdf[NUM_FREQS];
  Double_t rChiSquaresFullRange[NUM_FREQS];
  Int_t rNdfFullRange[NUM_FREQS];
  Double_t xHigh[NUM_FREQS];
  Double_t numOutliers[NUM_FREQS];

  Int_t nfreqs = NUM_FREQS;

  const int NUM_CHANS=500;
  const int CHAN_LOOPS=31; // total channels are something like 15500

  int maxchan, minchan;
  int loopchan;

  TFile *fout = new TFile(outFileName.c_str(), "recreate");
  TTree *rayleighTree = new TTree("rayleighTree", "Tree of Rayleigh fits of all channels");
  rayleighTree->Branch("run",           &run           , "run/I");
  rayleighTree->Branch("subrun",        &iev           , "subrun/I");
  rayleighTree->Branch("channel",       &channel       , "channel/I");
  rayleighTree->Branch("nfreqs",        &nfreqs        , "nfreqs/I");
  rayleighTree->Branch("powX",          powX           , "powX[nfreqs]/D");
  rayleighTree->Branch("rAmplitudes",   rAmplitudes    , "rAmplitudes[nfreqs]/D");
  rayleighTree->Branch("avgPowSpecY",   avgPowSpecY    , "avgPowSpecY[nfreqs]/D");
  rayleighTree->Branch("avgDiffPhases", avgDiffPhases  , "avgDiffPhases[nfreqs]/D"); 
  rayleighTree->Branch("rmsDiffPhases", rmsDiffPhases  , "rmsDiffPhases[nfreqs]/D"); 
  rayleighTree->Branch("rChiSquares",   rChiSquares    , "rChiSquares[nfreqs]/D");
  rayleighTree->Branch("rNdf",          rNdf           , "rNdf[nfreqs]/I");
 



  AveragePowerSpectrum *avgPowSpec = new AveragePowerSpectrum("avgPowSpec", "ProtoDUNE average power spectrum");   
    
  
  for (gallery::Event ev(filenames); !ev.atEnd(); ev.next()) {
    
    if(iev>=nevents) break;
    
    std::cout << "Event " << iev << std::endl;
    
    // Look at the digits                                                                                                                          
    auto& digits =
      *ev.getValidHandle<std::vector<raw::RawDigit>>(daq_tag);
    
    for(auto&& digit: digits){
      if(digit.Compression()!=0){
	std::cout << "Compression type " << digit.Compression() << std::endl;
      }
      // std::cout << iev << " " << digit.Channel() <<  " ";
      loopchan=digit.Channel();
      if ( loopchan != channel ) continue;
      
      
      mean=0;
      countsamples=0;
      for(auto&& sample: digit.ADCs()){
	// std::cout << sample << " ";
	samples[countsamples]=(sample);
	mean+=sample;
	countsamples++;
	if (countsamples>=NUM_SAMPLES) break;
      }
      mean/=(countsamples*1.);      
      
      for (int i=0; i<NUM_SAMPLES; i++) samples[i]-=mean;
      TGraph *gtemp = new TGraph(NUM_SAMPLES, x, samples);
      
      //      cout << "Event and channel " << iev << " " << channel << " " << mean << " "<< countsamples << " " << NUM_SAMPLES << endl;
      
      avgPowSpec->add(gtemp);
      
      if (doneGraph==false){
	sampleGraph = new TGraph(NUM_SAMPLES, x, samples);
	sampleFFT = FFTtools::makeRawPowerSpectrum(sampleGraph);

	doneGraph=true;
      }
      
      delete gtemp;
      
    } // end loop over digits (=?channels)                                                                                                        
    
    ++iev;
    
  } // end loop over events                                                                                                                        
  

  
  //  avgPowSpec->rebinAllRayleighHistograms(2); // emperically determined... for now.                                                                      
  avgPowSpec->fitAllRayleighHistograms();

  memset(powX, 0, sizeof(powX));
  memset(rAmplitudes, 0, sizeof(rAmplitudes));
  memset(rChiSquares, 0, sizeof(rChiSquares));
  memset(rNdf, 0, sizeof(rNdf));
  memset(rChiSquaresFullRange, 0, sizeof(rChiSquaresFullRange));
  memset(rNdfFullRange, 0, sizeof(rNdfFullRange));

  for (int ifreq=0; ifreq<NUM_FREQS; ifreq++){

    powX[ifreq] = ifreq*deltaf;

    avgPowSpecY[ifreq] = avgPowSpec->summedPowSpec[ifreq]/(avgPowSpec->count*1.);
    avgDiffPhases[ifreq] = avgPowSpec->summedDifferentialPhases[ifreq]/(avgPowSpec->count*1.);
    rmsDiffPhases[ifreq] = TMath::Sqrt(avgPowSpec->summedDifferentialPhasesSq[ifreq])/(avgPowSpec->count*1.);
    rAmplitudes[ifreq] = avgPowSpec->rayleighAmplitudes[ifreq];
    rChiSquares[ifreq] = avgPowSpec->rayleighFitChiSquares[ifreq];
    rNdf[ifreq] = avgPowSpec->rayleighNdf[ifreq];
    rChiSquaresOverNdf[ifreq] = rChiSquares[ifreq]*1./rNdf[ifreq];
    rChiSquaresFullRange[ifreq] = avgPowSpec->rayleighFitChiSquaresFullRange[ifreq];
    rNdfFullRange[ifreq] = avgPowSpec->rayleighNdfFullRange[ifreq];
     
    xHigh[ifreq] = avgPowSpec->xHigh[ifreq];
    numOutliers[ifreq] = avgPowSpec->numOutliers[ifreq];
    //      cout << " What number is this " << powX[ifreq]  << " " << summedPowSpec[ifreq]  << " " << rAmplitudes[ifreq] << endl;
  } 
      
  rayleighTree->Fill();
     
  fout->cd();
  cout << "Autosaving ..." << endl;
  rayleighTree->AutoSave();  

  sampleGraph->Write("sampleGraph");
  sampleFFT->Write("sampleFFT");

  TGraph *gAvgPower = new TGraph(NUM_FREQS, powX, avgPowSpecY);
  gAvgPower->SetTitle("Average power spectrum;Frequency [MHz];Amplitude");
  gAvgPower->Write("gAvgPowSpectrum");

  TGraph *gRayleighAmpl = new TGraph(NUM_FREQS, powX, rAmplitudes);
  gRayleighAmpl->SetTitle("Rayleigh Amplitude;Frequency [MHz];Amplitude [ADC/MHz]");
  gRayleighAmpl->Write("gRayleighAmpl");

  TGraph *gChiSquare = new TGraph(NUM_FREQS, powX, rChiSquaresOverNdf);
  gChiSquare->SetTitle("Rayleigh fit reduced chi square;Frequency [MHz];Chi square / NDF");
  gChiSquare->Write("gChiSquare");

  fout->cd();
  TDirectory *adir = gDirectory->mkdir("RayleighDistributions");
  adir->cd();

  for (int i=0; i<NUM_FREQS; i++){
    avgPowSpec->getRayleighHistogram(i)->Write();
    avgPowSpec->getRayleighHistogramFit(i)->Write();
    // avgPowSpec->getDiffPhaseHistogram(i)->Write();
  }


  // fout->Close();
  // delete avgPowSpec;  
  // cout << "Done and closed" << endl;

  // delete rayleighTree;
  // delete fout;
  // delete sampleGraph;

  cout << "Deleting everything now just need to return" << endl;

  return 0;

}

