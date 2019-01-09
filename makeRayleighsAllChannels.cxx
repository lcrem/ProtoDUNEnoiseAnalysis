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
#include <omp.h>

using namespace art;
using namespace std;
using namespace std::chrono;
//using namespace Acclaim;

void getChanInfo(UInt_t offlineChannel, UShort_t &apa, UShort_t &femb, UShort_t &fembch, UShort_t &plane, UShort_t &isOutside);

const int NUM_CHANS=200;
const int CHAN_LOOPS=11; // 77; // total offlineChannels is 2560*6=15360

void doEvent(int iev, gallery::Event &ev, InputTag const daq_tag, double samples[NUM_SAMPLES], AveragePowerSpectrum *avgPowSpec[NUM_CHANS], int iloop, double x[NUM_SAMPLES]);

int main(int argc, char *argv[]){

  UInt_t run;
  UInt_t subrun;
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

  UInt_t offlineChannel;

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


  int maxchan, minchan;
  int loopchan;
  UShort_t apa, apaChannel, femb, fembChannel;
  UShort_t isGood, isOutside, plane;

  TFile *fout = new TFile(outFileName.c_str(), "recreate");
  TTree *rayleighTree = new TTree("rayleighTree", "Tree of Rayleigh fits of all channels");
  rayleighTree->Branch("run",           &run           , "run/i");
  rayleighTree->Branch("subrun",        &subrun        , "subrun/i");
  rayleighTree->Branch("offlineChannel",&offlineChannel, "offlineChannel/i");
  rayleighTree->Branch("apa",           &apa           , "apa/s");
  rayleighTree->Branch("apaChannel",    &apaChannel    , "apaChannel/s");
  rayleighTree->Branch("femb",          &femb          , "femb/s");
  rayleighTree->Branch("fembChannel",   &fembChannel   , "fembChannel/s");
  rayleighTree->Branch("isGood",        &isGood        , "isGood/s");
  rayleighTree->Branch("isOutside",     &isOutside     , "isOutside/s");
  rayleighTree->Branch("plane",         &plane         , "plane/s");
  rayleighTree->Branch("nfreqs",        &nfreqs        , "nfreqs/I");
  rayleighTree->Branch("powX",          powX           , "powX[nfreqs]/D");
  rayleighTree->Branch("rAmplitudes",   rAmplitudes    , "rAmplitudes[nfreqs]/D");
  rayleighTree->Branch("avgPowSpecY",   avgPowSpecY    , "avgPowSpecY[nfreqs]/D");
  rayleighTree->Branch("avgDiffPhases", avgDiffPhases  , "avgDiffPhases[nfreqs]/D");
  rayleighTree->Branch("rmsDiffPhases", rmsDiffPhases  , "rmsDiffPhases[nfreqs]/D");
  rayleighTree->Branch("rChiSquares",   rChiSquares    , "rChiSquares[nfreqs]/D");
  rayleighTree->Branch("rNdf",          rNdf           , "rNdf[nfreqs]/I");
 

  // gallery::Event ev(filenames);
  // nevents = ev.numberOfEventsInFile() < nevents ? ev.numberOfEventsInFile() : nevents;


  struct timeval start, end;
  gettimeofday(&start, NULL);
 
  for (int iloop=0; iloop<CHAN_LOOPS; iloop++){

    // int this_thread = omp_get_thread_num(), num_threads = omp_get_num_threads();
    // cout << "Some info " << iloop << " " << this_thread << " " << num_threads << endl;
    
    minchan = iloop*NUM_CHANS;
    maxchan = (iloop+1)*NUM_CHANS;
    // loopchan=0;
    // mean=0;
    memset(samples, 0, sizeof(samples));


  AveragePowerSpectrum *avgPowSpec[NUM_CHANS];
  
  for (int ich=0; ich<NUM_CHANS; ich++) avgPowSpec[ich] = new AveragePowerSpectrum(Form("avgPowSpec_%d", ich), "ProtoDUNE average power spectrum");   
  
  // // run first event on one thread and then parallelise all other events
  // doEvent(0, ev, daq_tag, samples, avgPowSpec, iloop, x);
  // // #pragma omp parallel for  schedule(dynamic, 1)  private(samples)
  // for (iev=1;iev<nevents; ++iev){
  //   //      int this_thread = omp_get_thread_num();
  //   doEvent(iev, ev, daq_tag, samples, avgPowSpec, iloop, x);    
  // } // end loop over events                                                                                              
  
  iev=0;

  for (gallery::Event ev(filenames); !ev.atEnd(); ev.next()) {
    
    if(iev>=nevents){
      iev++;
      break;
    }
    std::cout << "Event " << iev << " in loop " << iloop << std::endl;
    
    // Look at the digits                                                                                                                          
    auto& digits =
      *ev.getValidHandle<std::vector<raw::RawDigit>>(daq_tag);
    
    for(auto&& digit: digits){
      if(digit.Compression()!=0){
	std::cout << "Compression type " << digit.Compression() << std::endl;
      }
      // std::cout << iev << " " << digit.Channel() <<  " ";
      // offlineChannel=digit.Channel();
      loopchan = digit.Channel() - (iloop*NUM_CHANS);
      
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
  
  cout << "Ended looping over events?" << endl;
  
  
  Double_t meanChi=0;
  // #pragma omp parallel for  schedule(dynamic, 1)  private(meanChi, offlineChannel, apa, apaChannel, femb, fembChannel, powX, rAmplitudes, rChiSquares, rNdf, isGood)
    for (int ich=0; ich<NUM_CHANS; ich++){

      //      avgPowSpec[ich]->rebinAllRayleighHistograms(2); // emperically determined... for now. 

      avgPowSpec[ich]->fitAllRayleighHistograms();

      // memset(powX, 0, sizeof(powX));
      // memset(rAmplitudes, 0, sizeof(rAmplitudes));
      // memset(rChiSquares, 0, sizeof(rChiSquares));
      // memset(rNdf, 0, sizeof(rNdf));
      // memset(rChiSquaresFullRange, 0, sizeof(rChiSquaresFullRange));
      // memset(rNdfFullRange, 0, sizeof(rNdfFullRange));

      offlineChannel = ich + iloop*NUM_CHANS;
      apa = (offlineChannel/2560);
      if (apa>6) continue;
      getChanInfo(offlineChannel, apaChannel, femb, fembChannel, plane, isOutside);
 
      for (int ifreq=0; ifreq<NUM_FREQS; ifreq++){

    	powX[ifreq] = ifreq*deltaf;

    	avgPowSpecY[ifreq] = avgPowSpec[ich]->summedPowSpec[ifreq]/(avgPowSpec[ich]->count*1.);
    	avgDiffPhases[ifreq] = avgPowSpec[ich]->summedDifferentialPhases[ifreq]/(avgPowSpec[ich]->count*1.);
	rmsDiffPhases[ifreq] = TMath::Sqrt(avgPowSpec[ich]->summedDifferentialPhasesSq[ifreq])/(avgPowSpec[ich]->count*1.);
    	rAmplitudes[ifreq] = avgPowSpec[ich]->rayleighAmplitudes[ifreq];
    	rChiSquares[ifreq] = avgPowSpec[ich]->rayleighFitChiSquares[ifreq];
    	rNdf[ifreq] = avgPowSpec[ich]->rayleighNdf[ifreq];
    	meanChi += rChiSquares[ifreq]*1.0/rNdf[ifreq];
    	// rChiSquaresOverNdf[ifreq] = rChiSquares[ifreq]/rNdf[ifreq];
    	// rChiSquaresFullRange[ifreq] = avgPowSpec[ich]->rayleighFitChiSquaresFullRange[ifreq];
    	// rNdfFullRange[ifreq] = avgPowSpec[ich]->rayleighNdfFullRange[ifreq];
    
    	// xHigh[ifreq] = avgPowSpec[ich]->xHigh[ifreq];
    	// numOutliers[ifreq] = avgPowSpec[ich]->numOutliers[ifreq];
    	//      cout << " What number is this " << powX[ifreq]  << " " << summedPowSpec[ifreq]  << " " << rAmplitudes[ifreq] << endl;
      } 
  
      if (meanChi/NUM_FREQS>10.) isGood=0;
      else isGood=1;

      rayleighTree->Fill();

      cout << "Done channel " << offlineChannel << " and loop is ( " << iloop << " ) and thread " << omp_get_thread_num()  << endl;

      delete avgPowSpec[ich];
    }

       
    fout->cd();
    cout << "Autosaving ..." << endl;
    rayleighTree->AutoSave();  


  }
  fout->Close();
  


  gettimeofday(&end, NULL);

  Int_t delta = (end.tv_sec  - start.tv_sec) ;

  cout << "Done and closed in " << delta << " seconds " << endl;

  return 0;

}


void getChanInfo(UInt_t offlineChannel, UShort_t &apaChannel, UShort_t &femb, UShort_t &fembch, UShort_t &plane, UShort_t &isOutside){

  UShort_t apa = offlineChannel/2560;

  apaChannel = offlineChannel%2560;
  
  isOutside=0; // induction plane

  // FEMBChannel ranges (u, v, x)
  //   X01 400-439, 1599-1560, 2080-2127
  //   X10 760-799, 1239-1200, 2512-2559
  //   X11 0-39, 1199-1160, 2079-2032
  //   X20 360-399, 839-800, 1647-1600

    if (apaChannel < 800){ // u plane
      plane = 0;
      if (apaChannel < 400){
	femb   = apaChannel/40 + 10;
	fembch = apaChannel%40;
      }else{ 
	femb   = (apaChannel-400)/40 ;
	fembch = (apaChannel-400)%40 ;
      }
    } else if (apaChannel < 1600){ // v plane
      plane = 1;
      if (apaChannel < 1200){
	femb   = 19 - (apaChannel-800)/40;
	fembch = 40 - (apaChannel-800)%40;
      }else{
	femb   = 9 - (apaChannel-1200)/40;
	fembch = 40 - (apaChannel-1200)%40;
      }    
    } else if (apaChannel < 2560 ) { // x plane
      plane = 2;
      if (apaChannel < 2080){
	femb   = 19 - (apaChannel-1600)/48;
	fembch = 48 - (apaChannel-1600)%48;
      } else {
	femb   = (apaChannel-2080)/48 ;
	fembch = (apaChannel-2080)%48 ;
      }

      // For even-numbered APAs, the wall-facing collection wires are 0-480; for odd-numbered APAs, 480-960
      if (apa%2==0 && (apaChannel-1600)<480 ) isOutside=1;
      else if (apa%2!=0 && (apaChannel-1600)>=480 ) isOutside=1;

    }



}



void doEvent(int iev, gallery::Event &ev, InputTag const daq_tag, double samples[NUM_SAMPLES], AveragePowerSpectrum *avgPowSpec[NUM_CHANS], int iloop, double x[NUM_SAMPLES]){

      ev.goToEntry(iev);
            
      std::cout << "Event " << iev << " in loop " << iloop << " and thread " << omp_get_thread_num() << std::endl;

      std::vector<raw::RawDigit> digits;
#pragma omp critical (createdigits)
      digits = *ev.getValidHandle<std::vector<raw::RawDigit>>(daq_tag);
      // auto const &digits =
      // 	*ev.getValidHandle<std::vector<raw::RawDigit>>(daq_tag);
	
      for(auto&& digit: digits){
	if(digit.Compression()!=0){
	  std::cout << "Compression type " << digit.Compression() << std::endl;
	}
	//	std::cout << iev << " " << digit.Channel() <<  " ";
	//	offlineChannel=digit.Channel();
	int loopchan = digit.Channel() - (iloop*NUM_CHANS);

	if ( loopchan<0 || loopchan >= NUM_CHANS ) continue;

	double mean=0;
	int countsamples=0;
	for(auto&& sample: digit.ADCs()){
	  // std::cout << sample << " ";
	  samples[countsamples]=(sample);
	  mean+=sample;
	  countsamples++;
	  if (countsamples>=NUM_SAMPLES) break;
	}
	mean/=(countsamples*1.);      
      

	for (int i=0; i<NUM_SAMPLES; i++) samples[i]-=mean;
	
	//      cout << "Event and offlineChannel " << iev << " " << offlineChannel << " " << mean << " "<< countsamples << " " << NUM_SAMPLES << endl;
	TGraph *gtemp = new TGraph(NUM_SAMPLES, x, samples);
	
	//#pragma omp critical (addtoavg)
	  avgPowSpec[loopchan]->add(gtemp);

	  // if (doneGraph==false){
	  //   sampleGraph = new TGraph(NUM_SAMPLES, x, samples);
	  //   doneGraph=true;
	  // }
	  
	  delete gtemp;
		
      } // end loop over digits (=?offlineChannels)  
      
}

