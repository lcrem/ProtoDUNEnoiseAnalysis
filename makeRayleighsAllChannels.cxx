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
    std::cerr << "Usage 1: " << argv[0] << " [irun] [subrun] (inlist) (outfile) (maxEvents)" << std::endl;
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

  int channel;

  double x[NUM_SAMPLES];
  
  double dt = NOMINAL_SAMPLING_DELTAT;
  double deltaf = 2./NUM_SAMPLES;

  for (int i=0; i<NUM_SAMPLES; i++) x[i] = i*dt;

  TGraph *sampleGraph;
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
  Double_t summedPowSpec[NUM_FREQS];
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
  rayleighTree->Branch("summedPowSpec", summedPowSpec  , "summedPowSpec[nfreqs]/D");
  rayleighTree->Branch("rChiSquares",   rChiSquares    , "rChiSquares[nfreqs]/D");
  rayleighTree->Branch("rNdf",          rNdf           , "rNdf[nfreqs]/I");
 



  for (int iloop=0; iloop<CHAN_LOOPS; iloop++){

    minchan = iloop*NUM_CHANS;
    maxchan = (iloop+1)*NUM_CHANS;

    AveragePowerSpectrum *avgPowSpec[NUM_CHANS];
    for (int ich=0; ich<NUM_CHANS; ich++) avgPowSpec[ich] = new AveragePowerSpectrum(Form("avgPowSpec_%d", ich), "ProtoDUNE average power spectrum");   
    
    
    iev=0;
    for (gallery::Event ev(filenames); !ev.atEnd(); ev.next()) {
      
      if(iev>=nevents) break;
      
      std::cout << "Event " << iev << " in loop " << iloop << std::endl;
            
      // Look at the digits                                                                                                                          
      auto& digits =
        *ev.getValidHandle<std::vector<raw::RawDigit>>(daq_tag);
      
      for(auto&& digit: digits){
	if(digit.Compression()!=0){
	  std::cout << "Compression type " << digit.Compression() << std::endl;
	}
	// std::cout << iev << " " << digit.Channel() <<  " ";
	channel=digit.Channel();
	loopchan = channel - (iloop*NUM_CHANS);

	if ( loopchan<0 || loopchan >= NUM_CHANS ) continue;


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

	avgPowSpec[loopchan]->add(gtemp);
      
	if (doneGraph==false){
	  sampleGraph = new TGraph(NUM_SAMPLES, x, samples);
	  doneGraph=true;
	}
      
	delete gtemp;

      } // end loop over digits (=?channels)                                                                                                        
    
      ++iev;

    } // end loop over events                                                                                                                        


    for (int ich=0; ich<NUM_CHANS; ich++){

      avgPowSpec[ich]->rebinAllRayleighHistograms(2); // emperically determined... for now.                                                                      
      avgPowSpec[ich]->fitAllRayleighHistograms();

      memset(powX, 0, sizeof(powX));
      memset(rAmplitudes, 0, sizeof(rAmplitudes));
      memset(rChiSquares, 0, sizeof(rChiSquares));
      memset(rNdf, 0, sizeof(rNdf));
      memset(rChiSquaresFullRange, 0, sizeof(rChiSquaresFullRange));
      memset(rNdfFullRange, 0, sizeof(rNdfFullRange));

      channel = ich + iloop*NUM_CHANS;


      for (int ifreq=0; ifreq<NUM_FREQS; ifreq++){

	powX[ifreq] = ifreq*deltaf;

	summedPowSpec[ifreq] = avgPowSpec[ich]->summedPowSpec[ifreq];
	rAmplitudes[ifreq] = avgPowSpec[ich]->rayleighAmplitudes[ifreq];
	rChiSquares[ifreq] = avgPowSpec[ich]->rayleighFitChiSquares[ifreq];
	rNdf[ifreq] = avgPowSpec[ich]->rayleighNdf[ifreq];
	rChiSquaresOverNdf[ifreq] = rChiSquares[ifreq]/rNdf[ifreq];
	rChiSquaresFullRange[ifreq] = avgPowSpec[ich]->rayleighFitChiSquaresFullRange[ifreq];
	rNdfFullRange[ifreq] = avgPowSpec[ich]->rayleighNdfFullRange[ifreq];
     
	xHigh[ifreq] = avgPowSpec[ich]->xHigh[ifreq];
	numOutliers[ifreq] = avgPowSpec[ich]->numOutliers[ifreq];
	//      cout << " What number is this " << powX[ifreq]  << " " << summedPowSpec[ifreq]  << " " << rAmplitudes[ifreq] << endl;
      } 
  
      rayleighTree->Fill();

      cout << "Done channel " << channel << " and loop is ( " << iloop << " ) " << endl;

      delete avgPowSpec[ich];
    }
    fout->cd();
    cout << "Autosaving ..." << endl;
    rayleighTree->AutoSave();  
  }


  fout->Close();
  
  cout << "Done and closed" << endl;

  delete rayleighTree;
  delete fout;
  delete sampleGraph;

  cout << "Deleting everything now just need to return" << endl;

  return 0;

}

