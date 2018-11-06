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
  string inFileName;
  string outFileName;
  Int_t nevents=1e9;
  
  if((argc<2)){
    std::cerr << "Usage 1: " << argv[0] << " [irun] [isubrun] (infile) (outfile) (maxEvents)" << std::endl;
    return 1;
  } else {
    run = atoi(argv[1]);
    subrun = atoi(argv[2]);
    if (argc>3) inFileName += argv[3];
    else  inFileName +=  "../filenameList.txt";
    if (argc>4) outFileName += argv[4];
    else outFileName = "tempRayleighTree.root";
    if (argc>5)  nevents = atoi(argv[5]);

  }

  cout << " Processing " << inFileName << endl;

  int offlineChannel;

  double x[10000];
  
  double dt = NOMINAL_SAMPLING_DELTAT;
  double deltaf = 2./NUM_SAMPLES;

  for (int i=0; i<NUM_SAMPLES; i++) x[i] = i*dt;

  TGraph *sampleGraph;
  bool doneGraph=false;

  InputTag const daq_tag{ "tpcrawdecoder", "daq", "DecoderandReco" }; 

  // Create a vector of length 1, containing the given filename.                                                                                   
  double samples[10000];

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


  TFile *fout = new TFile(outFileName.c_str(), "recreate");
  TTree *rayleighTree = new TTree("rayleighTree", "Tree of Rayleigh fits of all channels");
  
  Double_t powX[NUM_FREQS];
  Double_t rAmplitudes[NUM_FREQS];
  Double_t summedPowSpec[NUM_FREQS];
  Double_t rChiSquares[NUM_FREQS];
  Double_t rChiSquaresOverNdf[NUM_FREQS];
  Int_t rNdf[NUM_FREQS];
  Double_t rChiSquaresFullRange[NUM_FREQS];
  Int_t rNdfFullRange[NUM_FREQS];
  Double_t xHigh[NUM_FREQS];
  Double_t numOutliers[NUM_FREQS];

  Int_t nfreqs = NUM_FREQS;

  rayleighTree->Branch("run",           &run           , "run/I");
  rayleighTree->Branch("subrun",        &subrun        , "subrun/I");
  rayleighTree->Branch("iev",           &iev           , "iev/I");
  rayleighTree->Branch("nfreqs",        &nfreqs        , "nfreqs/I");
  rayleighTree->Branch("powX",          powX           , "powX[nfreqs]/D");
  rayleighTree->Branch("rAmplitudes",   rAmplitudes    , "rAmplitudes[nfreqs]/D");
  rayleighTree->Branch("summedPowSpec", summedPowSpec  , "summedPowSpec[nfreqs]/D");
  rayleighTree->Branch("rChiSquares",   rChiSquares    , "rChiSquares[nfreqs]/D");
  rayleighTree->Branch("rNdf",          rNdf           , "rNdf[nfreqs]/I");
 

  for (gallery::Event ev(filenames); !ev.atEnd(); ev.next()) {
    
    if(iev>=nevents) break;

    std::cout << "Event " << iev << std::endl;
    
    AveragePowerSpectrum *avgPowSpec = new AveragePowerSpectrum("avgPowSpec", "ProtoDUNE average power spectrum");   
    // Look at the digits                                                                                                                          
    auto& digits =
        *ev.getValidHandle<std::vector<raw::RawDigit>>(daq_tag);

    for(auto&& digit: digits){
      if(digit.Compression()!=0){
        std::cout << "Compression type " << digit.Compression() << std::endl;
      }
      // std::cout << iev << " " << digit.Channel() <<  " ";
      offlineChannel=digit.Channel();
      mean=0;
      countsamples=0;
      for(auto&& sample: digit.ADCs()){
	// std::cout << sample << " ";
	samples[countsamples]=(sample);
	mean+=sample;
	countsamples++;
      }
      mean/=(countsamples*1.);      
      
      for (int i=0; i<NUM_SAMPLES; i++) samples[i]-=mean;
      TGraph *gtemp = new TGraph(NUM_SAMPLES, x, samples);
      
      avgPowSpec->add(gtemp);
      //      cout << "Event and channel " << iev << " " << offlineChannel << " " << mean << " "<< countsamples << " " << NUM_SAMPLES << endl;
      
      if (doneGraph==false){
	sampleGraph = new TGraph(NUM_SAMPLES, x, samples);
	doneGraph=true;
      }
      
      delete gtemp;

    } // end loop over digits (=?channels)                                                                                                        
    

    
    avgPowSpec->rebinAllRayleighHistograms(2); // emperically determined... for now.                                                                      
    avgPowSpec->fitAllRayleighHistograms();

    memset(powX, 0, sizeof(powX));
    memset(rAmplitudes, 0, sizeof(rAmplitudes));
    memset(rChiSquares, 0, sizeof(rChiSquares));
    memset(rNdf, 0, sizeof(rNdf));
    memset(rChiSquaresFullRange, 0, sizeof(rChiSquaresFullRange));
    memset(rNdfFullRange, 0, sizeof(rNdfFullRange));


    for (int ifreq=0; ifreq<NUM_FREQS; ifreq++){

      powX[ifreq] = ifreq*deltaf;

      summedPowSpec[ifreq] = avgPowSpec->summedPowSpec[ifreq];
      rAmplitudes[ifreq] = avgPowSpec->rayleighAmplitudes[ifreq];
      rChiSquares[ifreq] = avgPowSpec->rayleighFitChiSquares[ifreq];
      rNdf[ifreq] = avgPowSpec->rayleighNdf[ifreq];
      rChiSquaresOverNdf[ifreq] = rChiSquares[ifreq]/rNdf[ifreq];
      rChiSquaresFullRange[ifreq] = avgPowSpec->rayleighFitChiSquaresFullRange[ifreq];
      rNdfFullRange[ifreq] = avgPowSpec->rayleighNdfFullRange[ifreq];
     
      xHigh[ifreq] = avgPowSpec->xHigh[ifreq];
      numOutliers[ifreq] = avgPowSpec->numOutliers[ifreq];
      //      cout << " What number is this " << powX[ifreq]  << " " << summedPowSpec[ifreq]  << " " << rAmplitudes[ifreq] << endl;
    } 
  


    rayleighTree->Fill();
    ++iev;


    //avgPowSpec->deleteRayleighDistributions();

    delete avgPowSpec;

  } // end loop over events                                                                                                                        

  // TDirectory *adir = gDirectory->mkdir("RayleighDistributions");
  // adir->cd();

  // for (int i=0; i<NUM_FREQS; i++){
  //   avgPowSpec->getRayleighHistogram(i)->Write();
  //   avgPowSpec->getRayleighHistogramFit(i)->Write();
  // }

  fout->cd();
  rayleighTree->Write();
  fout->Close();

  return 0;

}

