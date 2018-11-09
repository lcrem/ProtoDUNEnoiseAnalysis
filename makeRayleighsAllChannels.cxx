#include "AveragePowerSpectrum.h"
#include "canvas/Utilities/InputTag.h"
#include "gallery/Event.h"
#include "lardataobj/RawData/RawDigit.h"
#include "dune-raw-data/Services/ChannelMap/PdspChannelMapService.h"

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

int getFEMBFromOfflineChannel(int offlineChannel);
int getFEMBChannelFromOfflineChannel(int offlineChannel);

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

  int offlineChannel;

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

  const int NUM_CHANS=100;
  const int CHAN_LOOPS=154; // total offlineChannels is 2560*6=15360

  int maxchan, minchan;
  int loopchan;
  int apa, apaChannel, femb, fembChannel;

  TFile *fout = new TFile(outFileName.c_str(), "recreate");
  TTree *rayleighTree = new TTree("rayleighTree", "Tree of Rayleigh fits of all channels");
  rayleighTree->Branch("run",           &run           , "run/I");
  rayleighTree->Branch("subrun",        &iev           , "subrun/I");
  rayleighTree->Branch("offlineChannel",&offlineChannel, "offlineChannel/I");
  rayleighTree->Branch("apa",           &apa           , "apa/I");
  rayleighTree->Branch("apaChannel",    &apaChannel    , "apaChannel/I");
  rayleighTree->Branch("femb",          &femb          , "femb/I");
  rayleighTree->Branch("fembChannel",   &fembChannel   , "fembChannel/I");
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
	offlineChannel=digit.Channel();
	loopchan = offlineChannel - (iloop*NUM_CHANS);

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
      
	//      cout << "Event and offlineChannel " << iev << " " << offlineChannel << " " << mean << " "<< countsamples << " " << NUM_SAMPLES << endl;

	avgPowSpec[loopchan]->add(gtemp);
      
	if (doneGraph==false){
	  sampleGraph = new TGraph(NUM_SAMPLES, x, samples);
	  doneGraph=true;
	}
      
	delete gtemp;

      } // end loop over digits (=?offlineChannels)                                                                                                        
    
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

      offlineChannel = ich + iloop*NUM_CHANS;
      apa = (offlineChannel/2560);
      apaChannel = (offlineChannel%2560);
      femb = getFEMBFromOfflineChannel(offlineChannel);
      fembChannel = getFEMBChannelFromOfflineChannel(offlineChannel);

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

      cout << "Done channel " << offlineChannel << " and loop is ( " << iloop << " ) " << endl;

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

int getFEMBFromOfflineChannel(int offlineChannel){

  int femb=-1;

  int apaChannel = offlineChannel%2560;

  // FEMBChannel ranges (u, v, x)
  //   X01 400-439, 1599-1560, 2080-2127
  //   X10 760-799, 1239-1200, 2512-2559
  //   X11 0-39, 1199-1160, 2079-2032
  //   X20 360-399, 839-800, 1647-1600

    if (apaChannel < 800){ // u plane
      if (apaChannel < 400)
	femb = apaChannel/40 + 11;
      else 
	femb = (apaChannel-400)/40 + 1;
      
    } else if (apaChannel < 1600){ // v plane
      if (apaChannel < 1200)
	femb = 20 - (apaChannel-800)/40;
      else
	femb = 10 - (apaChannel-1200)/40;
    
    } else if (apaChannel < 2560 ) { // x plane
      if (apaChannel < 2080)
	femb = 20 - (apaChannel-1600)/48;
      else
	femb = (apaChannel-2080)/48  + 1 ;
    }

  return femb;
}


int getFEMBChannelFromOfflineChannel(int offlineChannel){

  int apaChannel = offlineChannel%2560;
  
  int fembch = -1;

  // FEMBChannel ranges (u, v, x)
  //   X01 400-439, 1599-1560, 2080-2127
  //   X10 760-799, 1239-1200, 2512-2559
  //   X11 0-39, 1199-1160, 2079-2032
  //   X20 360-399, 839-800, 1647-1600

    if (apaChannel < 800){ // u plane
      if (apaChannel < 400)
	fembch = apaChannel%40;
      else 
	fembch = (apaChannel-400)%40 ;
      
    } else if (apaChannel < 1600){ // v plane
      if (apaChannel < 1200)
	fembch = 40 - (apaChannel-800)%40;
      else
	fembch = 40 - (apaChannel-1200)%40;
    
    } else if (apaChannel < 2560 ) { // x plane
      if (apaChannel < 2080)
	fembch = 48 - (apaChannel-1600)%48;
      else
	fembch = (apaChannel-2080)%48 ;
    }

  return fembch;

}
